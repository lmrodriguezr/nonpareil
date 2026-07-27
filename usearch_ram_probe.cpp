// usearch_ram_probe - Standalone experiment: measure how usearch's
// -makeudb_usearch peak RSS actually scales with input size, since a
// 2-point linear extrapolation (from ~2-8k read calibration builds) was
// observed to underestimate real production memory by >4x on a shard with
// millions of reads. Not part of the nonpareil pipeline itself -- this is
// a diagnostic tool to characterize the real shape of the function before
// changing the production calibration/sharding logic.
//
// @license artistic 2.0
//
// Build: `make usearch_ram_probe` from the repo root.
//
// Usage -- scaling curve (peak RSS vs read count/bases, at increasing N,
// using the same -slots formula (2*N) the production code uses):
//   ./usearch_ram_probe scale <input.fastq[.gz]> <fasta|fastq> <out.csv> \
//     <N1> <N2> ...
//   e.g. ./usearch_ram_probe scale my.fastq.gz fastq scale.csv \
//     2000 8000 32000 128000 512000 2000000 6000000
//
// Usage -- slots sensitivity at a fixed read count (isolates whether
// -slots/hash-table sizing drives RSS independently of bases/reads):
//   ./usearch_ram_probe slots <input.fastq[.gz]> <fasta|fastq> <out.csv> \
//     <fixedN> <mult1> <mult2> ...
//   e.g. ./usearch_ram_probe slots my.fastq.gz fastq slots.csv \
//     200000 0.25 0.5 1 2 4 8
//
// Both modes stream results to <out.csv> after every build (so partial
// results survive if a later, larger build gets OOM-killed), and the
// scale sweep stops automatically once a build is killed by a signal
// (larger N would almost certainly fail too).

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <sys/wait.h>
#include <sys/resource.h>
#include <unistd.h>
#include <filesystem>

#include "enveomics/universal.h"
#include "enveomics/sequence.h"

#define LARGEST_PATH 4096

struct RunResult {
  bool   ran = false;
  int    exit_code = -1;
  bool   signaled = false;
  int    term_signal = 0;
  long   peak_rss_kb = -1;
  double wall_seconds = 0;
};

// Runs one child process to completion and returns *its own* peak RSS via
// wait4(), rather than the getrusage(RUSAGE_CHILDREN) baseline-delta trick
// the production calibration code uses -- that trick can misreport when
// builds aren't strictly increasing in memory use (RUSAGE_CHILDREN's
// maxrss is a running max across all children, so a delta reads ~0 if a
// later child happens to peak lower than an earlier one). wait4() gives an
// exact, unambiguous answer for this one process, with no ordering
// assumptions.
static RunResult run_and_measure(const std::vector<std::string> &args) {
  RunResult r;
  std::vector<char*> argv;
  for (auto &s : args) argv.push_back(const_cast<char*>(s.c_str()));
  argv.push_back(nullptr);

  auto t0 = std::chrono::steady_clock::now();
  pid_t pid = fork();
  if (pid < 0) error("fork() failed");
  if (pid == 0) {
    execvp(argv[0], argv.data());
    perror("execvp");
    _exit(127);
  }

  int status;
  struct rusage ru;
  pid_t w = wait4(pid, &status, 0, &ru);
  auto t1 = std::chrono::steady_clock::now();
  r.wall_seconds = std::chrono::duration<double>(t1 - t0).count();
  if (w == pid) {
    r.ran = true;
    if (WIFEXITED(status))   r.exit_code = WEXITSTATUS(status);
    if (WIFSIGNALED(status)) { r.signaled = true; r.term_signal = WTERMSIG(status); }
#ifdef __APPLE__
    r.peak_rss_kb = ru.ru_maxrss / 1024; // macOS reports bytes, Linux Kb
#else
    r.peak_rss_kb = ru.ru_maxrss;
#endif
  }
  return r;
}

static size_t sum_lengths(const std::vector<unsigned int> &lens, size_t n) {
  size_t sum = 0;
  for (size_t i = 0; i < n && i < lens.size(); i++) sum += lens[i];
  return sum;
}

int main(int argc, char **argv) {
  set_verbosity(4);
  if (argc < 6) {
    std::cerr
      << "Usage:\n"
      << "  " << argv[0] << " scale <input> <fasta|fastq> <out.csv> <N1> [N2 ...]\n"
      << "  " << argv[0] << " slots <input> <fasta|fastq> <out.csv> <fixedN> <mult1> [mult2 ...]\n";
    return 1;
  }

  std::string mode  = argv[1];
  char *inputFile   = argv[2];
  char *format      = argv[3];
  std::string outCsv = argv[4];

  std::vector<double> values;
  for (int i = 5; i < argc; i++) values.push_back(atof(argv[i]));

  // Builds the full index once, exactly like the production usearch
  // sharding path does -- reuses build_index (gz-aware) and
  // get_seq_lengths/write_seq_range_to_fasta, so "bases"/"reads" here mean
  // precisely what they mean in enveomics/nonpareil_mating.cpp.
  char *namFile, *seqFile;
  int largest_seq;
  double avg_seq;
  std::cerr << "Indexing " << inputFile << " ...\n";
  size_t total_seqs = build_index(
    inputFile, format, namFile, seqFile, largest_seq, avg_seq, 1, false
  );
  std::cerr << "Indexed " << total_seqs << " reads, avg length "
            << avg_seq << " bp, longest " << largest_seq << " bp\n";

  std::vector<unsigned int> lens = get_seq_lengths(
    seqFile, (unsigned int) total_seqs
  );
  std::filesystem::path tmp = tmp_dir();

  std::ofstream csv(outCsv);
  if (!csv.is_open()) error("Cannot open output CSV", outCsv.c_str());

  if (mode == "scale") {
    csv << "reads,bases,slots,peak_rss_kb,wall_seconds,exit_code,signaled,term_signal\n";
    for (double v : values) {
      size_t n = (size_t) v;
      if (n > total_seqs) {
        std::cerr << "Skipping N=" << n << " (only " << total_seqs
                   << " reads available)\n";
        continue;
      }
      size_t bases = sum_lengths(lens, n);
      size_t slots = n * 2; // matches production's default -slots formula

      char fastaPath[LARGEST_PATH], dbPath[LARGEST_PATH];
      snprintf(fastaPath, LARGEST_PATH, "%s/probe_%zu.fasta", tmp.c_str(), n);
      snprintf(dbPath,    LARGEST_PATH, "%s/probe_%zu.db",    tmp.c_str(), n);
      write_seq_range_to_fasta(seqFile, fastaPath, 1, n);

      std::cerr << "Building N=" << n << " reads (" << bases
                << " bases, slots=" << slots << ") ...\n";
      RunResult r = run_and_measure({
        "usearch", "-makeudb_usearch", fastaPath,
        "-output", dbPath, "-slots", std::to_string(slots)
      });

      csv << n << "," << bases << "," << slots << ","
          << r.peak_rss_kb << "," << r.wall_seconds << ","
          << r.exit_code << "," << (r.signaled ? 1 : 0) << ","
          << r.term_signal << "\n";
      csv.flush();

      std::cerr << "  -> peak_rss=" << r.peak_rss_kb << " Kb, wall="
                << r.wall_seconds << "s";
      if (r.signaled)
        std::cerr << "  [KILLED by signal " << r.term_signal << "]";
      std::cerr << "\n";

      remove(fastaPath);
      remove(dbPath);

      if (r.signaled) {
        std::cerr << "Stopping scale sweep: this build was killed (likely "
                   << "OOM) -- larger N would almost certainly fail too.\n";
        break;
      }
    }
  } else if (mode == "slots") {
    if (values.size() < 2)
      error("slots mode needs <fixedN> and at least one multiplier");
    size_t n = (size_t) values[0];
    if (n > total_seqs) error("Requested N exceeds available reads");
    size_t bases = sum_lengths(lens, n);

    char fastaPath[LARGEST_PATH];
    snprintf(fastaPath, LARGEST_PATH, "%s/probe_slots_%zu.fasta", tmp.c_str(), n);
    write_seq_range_to_fasta(seqFile, fastaPath, 1, n);

    csv << "reads,bases,slots,slot_multiplier,peak_rss_kb,wall_seconds,"
        << "exit_code,signaled,term_signal\n";
    for (size_t i = 1; i < values.size(); i++) {
      double mult = values[i];
      size_t slots = (size_t)(n * 2 * mult);
      if (slots < 1) slots = 1;

      char dbPath[LARGEST_PATH];
      snprintf(
        dbPath, LARGEST_PATH, "%s/probe_slots_%zu_%zu.db",
        tmp.c_str(), n, slots
      );

      std::cerr << "Building N=" << n << " reads with slots=" << slots
                << " (mult=" << mult << ") ...\n";
      RunResult r = run_and_measure({
        "usearch", "-makeudb_usearch", fastaPath,
        "-output", dbPath, "-slots", std::to_string(slots)
      });

      csv << n << "," << bases << "," << slots << "," << mult << ","
          << r.peak_rss_kb << "," << r.wall_seconds << ","
          << r.exit_code << "," << (r.signaled ? 1 : 0) << ","
          << r.term_signal << "\n";
      csv.flush();

      std::cerr << "  -> peak_rss=" << r.peak_rss_kb << " Kb\n";
      remove(dbPath);
    }
    remove(fastaPath);
  } else {
    error("Unknown mode (expected 'scale' or 'slots')", mode.c_str());
  }

  std::cerr << "Done. Results in " << outCsv << "\n";
  return 0;
}
