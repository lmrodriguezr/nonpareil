// nonpareil_mating - Part of the nonpareil package
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @license artistic 2.0

#define _MULTI_THREADED
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <sstream>
#include <filesystem>
#include <vector>
#include <algorithm>

#include "universal.h"
#include "multinode.h"
#include "sequence.h"
#include "nonpareil_mating.h"

extern int processID;
extern int processes;

#define LARGEST_PATH 4096
#define USEARCH_RAM_SAFETY_FACTOR 0.75
// -maxrejects 0 (fully exhaustive) turned out to be the dominant cost of
// the usearch kernel: for the (common) case of reads with few/no true
// hits, it forces scanning the entire remaining shard just to confirm
// there's nothing left. Bounding it is a deliberate, calibrated trade-off,
// not just a performance knob -- the coverage-from-redundancy model
// (`Nonpareil.kappa_to_coverage`'s `common_factor` in
// utils/Nonpareil/R/Nonpareil.R, for kernel=="usearch") must be refit
// against whatever value is used here, since it's what "redundancy" means
// as far as that model is concerned. Do not change this without
// recalibrating that constant to match.
#define USEARCH_MAXREJECTS 100
#define USEARCH_CALIB_SMALL 8000u
#define USEARCH_CALIB_LARGE 32000u

namespace {

// A contiguous, RAM-bounded slice of the subject dataset (1-based `start`,
// matching the indexing convention of `get_seqs` elsewhere in this file).
struct shard_t {
  size_t       start;
  unsigned int count;
  int          max_len;
  size_t       total_bases;
};

// Runs one child process to completion via fork()+execvp(), redirecting its
// stdout/stderr to `log_file`, and returns its *exact* peak RSS (Kb) via
// wait4(). This deliberately does NOT use system() + a
// getrusage(RUSAGE_CHILDREN) delta (the calibration code's original
// approach): RUSAGE_CHILDREN's ru_maxrss is a running max across every
// child since process start, not a per-call reading, so a delta taken
// around a *second* (or later) build silently measures
// (that build's peak - the previous build's peak) instead of an absolute
// value. That's exactly what made a real production run's calibration
// undermeasure its slope by ~3x (confirmed empirically -- see
// usearch_ram_probe.cpp/the calibration experiment this was validated
// against). wait4() reports the exact rusage of the one child just
// waited on, with no ordering or baseline assumptions -- do not
// "simplify" this back to system()+getrusage(RUSAGE_CHILDREN).
long run_and_get_peak_rss_kb(
      const std::vector<std::string> &args, const char *log_file) {
  std::vector<char*> argv;
  for (const std::string &s : args) argv.push_back(const_cast<char*>(s.c_str()));
  argv.push_back(NULL);

  pid_t pid = fork();
  if (pid < 0) error("fork() failed while launching usearch");
  if (pid == 0) {
    int fd = open(log_file, O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd >= 0) {
      dup2(fd, STDOUT_FILENO);
      dup2(fd, STDERR_FILENO);
      close(fd);
    }
    execvp(argv[0], argv.data());
    _exit(127); // exec failed
  }

  int status;
  struct rusage ru;
  pid_t w = wait4(pid, &status, 0, &ru);
  if (w != pid) error("wait4() failed while waiting for usearch");
  if (WIFSIGNALED(status))
    error("usearch (calibration) was killed by signal", WTERMSIG(status));
  if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
    error(
      "usearch 'makeudb' (calibration) failed with return code",
      WIFEXITED(status) ? WEXITSTATUS(status) : -1
    );

#ifdef __APPLE__
  return ru.ru_maxrss / 1024; // macOS reports bytes, Linux reports Kb
#else
  return ru.ru_maxrss;
#endif
}

// Builds a usearch index from the first `n` reads of `file`, and returns
// the exact peak RSS (Kb) of that single build.
long usearch_calibration_point(char *file, unsigned int n, size_t &out_bases) {
  std::filesystem::path base = tmp_dir();
  char *fasta = new char[LARGEST_PATH],
       *db    = new char[LARGEST_PATH],
       *log   = new char[LARGEST_PATH],
       *cmd_display = new char[LARGEST_PATH];
  snprintf(fasta, LARGEST_PATH, "%s/usearch_calib_%u.fasta", base.c_str(), n);
  snprintf(db,    LARGEST_PATH, "%s/usearch_calib_%u.db",    base.c_str(), n);
  snprintf(log,   LARGEST_PATH, "%s/usearch_calib_%u.log",   base.c_str(), n);

  write_seq_range_to_fasta(file, fasta, 1, n);

  std::vector<unsigned int> lens = get_seq_lengths(fasta, n);
  out_bases = 0;
  for (size_t a = 0; a < lens.size(); a++) out_bases += lens[a];

  snprintf(
    cmd_display, LARGEST_PATH,
    "usearch -makeudb_usearch '%s' -output '%s' -slots %u (log: '%s')",
    fasta, db, n * 2, log
  );
  say("5ss$", "CMD (calibration): ", cmd_display);
  long peak = run_and_get_peak_rss_kb(
    {
      "usearch", "-makeudb_usearch", fasta, "-output", db,
      "-slots", std::to_string(n * 2)
    },
    log
  );

  remove(fasta);
  remove(db);
  remove(log);
  delete[] fasta;
  delete[] db;
  delete[] log;
  delete[] cmd_display;

  return peak;
}

// Fits `peak_kb ~= a + b * total_bases` from two calibration points built
// from the start of `file` (sizes bounded by USEARCH_CALIB_SMALL/_LARGE
// regardless of dataset size, so calibration itself stays cheap), then walks
// the per-read lengths of `file` to produce a list of RAM-bounded,
// contiguous shards. `ram_mb_per_rank` is assumed to already be adjusted for
// any co-located MPI ranks.
std::vector<shard_t> plan_shards(
      char *file, unsigned int total_seqs, double ram_mb_per_rank) {
  std::vector<shard_t> shards;
  std::vector<unsigned int> lengths = get_seq_lengths(file, total_seqs);

  double bases_budget = -1;
  unsigned int calib_a = std::min(USEARCH_CALIB_SMALL, total_seqs);
  unsigned int calib_b = std::min(USEARCH_CALIB_LARGE, total_seqs);

  if (calib_a > 0 && calib_b > calib_a) {
    size_t bases_a, bases_b;
    long rss_a = usearch_calibration_point(file, calib_a, bases_a);
    long rss_b = usearch_calibration_point(file, calib_b, bases_b);

    if (bases_b > bases_a) {
      double b_coef = (double)(rss_b - rss_a) / (double)(bases_b - bases_a);
      double a_coef = (double) rss_a - b_coef * (double) bases_a;
      say("5sfsfs$", "Calibration: a=", a_coef, " Kb, b=", b_coef, " Kb/base");
      if (b_coef > 0)
        bases_budget = (
          ram_mb_per_rank * 1024.0 * USEARCH_RAM_SAFETY_FACTOR - a_coef
        ) / b_coef;
    }
  }

  if (bases_budget <= 0) {
    // Degenerate calibration (dataset too small for two distinct points, or
    // a non-positive fitted slope): fall back to a single shard spanning
    // the whole dataset, equivalent to the un-sharded orientation.
    size_t total_bases = 0;
    int    max_len = 0;
    for (size_t a = 0; a < lengths.size(); a++) {
      total_bases += lengths[a];
      if ((int) lengths[a] > max_len) max_len = lengths[a];
    }
    shard_t shard = {1, total_seqs, max_len, total_bases};
    shards.push_back(shard);
    return shards;
  }

  say("3sfs$", "Shard budget: ", bases_budget, " bases");

  size_t cur_start = 1, cur_bases = 0;
  unsigned int cur_count = 0;
  int cur_max_len = 0;
  for (unsigned int i = 0; i < total_seqs; i++) {
    unsigned int len = lengths[i];
    if (cur_count > 0 && (double)(cur_bases + len) > bases_budget) {
      shard_t shard = {cur_start, cur_count, cur_max_len, cur_bases};
      shards.push_back(shard);
      cur_start += cur_count;
      cur_count = 0;
      cur_bases = 0;
      cur_max_len = 0;
    }
    cur_count++;
    cur_bases += len;
    if ((int) len > cur_max_len) cur_max_len = len;
  }
  if (cur_count > 0) {
    shard_t shard = {cur_start, cur_count, cur_max_len, cur_bases};
    shards.push_back(shard);
  }

  return shards;
}

// Single sequential pass over `file`: writes each shard owned by this rank
// (`shard index % processes == processID`) to its own FASTA file at
// "<tmp_dir>/usearch_shard_<start>.fasta", skipping (not buffering) reads
// belonging to shards owned by other ranks. This avoids re-scanning the
// whole dataset once per shard, which would be O(shards x total_seqs).
void write_owned_shards_to_fasta(
      char *file, const std::vector<shard_t> &shards) {
  ifstream filein;
  ofstream fileout;
  string   entry, header;
  size_t   i = 0, shard_i = 0;
  bool     out_open = false;
  std::filesystem::path base = tmp_dir();

  filein.open(file, ios::in);
  if (!filein.is_open()) error("Impossible to open the input file", file);

  while (filein.good()) {
    string line;
    getline(filein, line);
    if ((line.size() > 0 && line[0] == '>') || !filein.good()) {
      if (entry.size() > 0) {
        i++;
        while (
              shard_i < shards.size() &&
              i > shards[shard_i].start + shards[shard_i].count - 1) {
          if (out_open) { fileout.close(); out_open = false; }
          shard_i++;
        }
        if (shard_i < shards.size() && i >= shards[shard_i].start &&
              ((int)(shard_i % processes) == processID)) {
          if (!out_open) {
            char outfile[LARGEST_PATH];
            snprintf(
              outfile, LARGEST_PATH, "%s/usearch_shard_%zu.fasta",
              base.c_str(), shards[shard_i].start
            );
            fileout.open(outfile, ios::out);
            if (!fileout.is_open())
              error("Impossible to open the output file", outfile);
            out_open = true;
          }
          fileout << header << "\n" << entry << "\n";
          if (fileout.fail())
            error("Write to shard fasta failed", (int) shard_i);
        }
      }
      header = line;
      entry = (string)"";
    } else {
      entry.append(line);
    }
  }
  if (out_open) fileout.close();
  filein.close();
}

// Builds the shard's usearch index, queries it with `sampleFile` (the
// subsample) as query -- the natural orientation, since the shard is now
// small enough to be the DB -- and folds matching hits into `result`.
// Deletes all of the shard's temporary files before returning.
void process_one_shard(
      int *&result, char *sampleFile, const shard_t &shard, int threads,
      size_t qry_seqs, matepar_t matepar) {
  std::filesystem::path base = tmp_dir();
  char fasta[LARGEST_PATH], db[LARGEST_PATH], out[LARGEST_PATH],
       log[LARGEST_PATH], cmd1[LARGEST_PATH], cmd2[LARGEST_PATH];
  snprintf(
    fasta, LARGEST_PATH, "%s/usearch_shard_%zu.fasta",
    base.c_str(), shard.start
  );
  snprintf(
    db, LARGEST_PATH, "%s/usearch_shard_%zu.db", base.c_str(), shard.start
  );
  snprintf(
    out, LARGEST_PATH, "%s/usearch_shard_%zu.out", base.c_str(), shard.start
  );
  snprintf(
    log, LARGEST_PATH, "%s/usearch_shard_%zu.log", base.c_str(), shard.start
  );

  // Index the shard's USearch DB
  size_t slots = matepar.hashsize;
  if (slots == 0) slots = (size_t)(shard.count * 2);
  snprintf(
    cmd1, LARGEST_PATH,
    "usearch -makeudb_usearch '%s' -output '%s' -slots %zu > '%s' 2>&1",
    fasta, db, slots, log
  );
  say("8ss$", "CMD: ", cmd1);
  int ret1 = system(cmd1);
  if (ret1 != 0) error("usearch 'makeudb' failed with return code", ret1);

  // Query the subsample against the shard
  snprintf(
    cmd2, LARGEST_PATH,
    "usearch -usearch_local '%s' -db '%s' -userout '%s' -threads '%d' \
      -evalue 0.00001 -id 0.9 -userfields '%s' -strand both \
      -maxaccepts 0 -maxrejects %d \
      >> '%s' 2>&1",
    sampleFile, db, out, threads, "query+target+qcov+tcov",
    USEARCH_MAXREJECTS, log
  );
  say("8ss$", "CMD: ", cmd2);
  int ret2 = system(cmd2);
  if (ret2 != 0) error("usearch 'local' failed with return code", ret2);

  // Parse the output.
  ifstream filein;
  filein.open(out, ios::in);
  if (!filein.is_open()) error("Impossible to open the input file", out);
  while (filein.good()) {
    string line;
    getline(filein, line);
    if (line.size() == 0) continue;

    std::vector<string> fields;
    string token;
    stringstream ss(line);
    while (getline(ss, token, '\t')) fields.push_back(token);
    if (fields.size() < 4) continue; // not enough columns

    try {
      int tid = stoi(fields[0]); // <- This is the subsample "query"
      double qcov = stod(fields[2]);
      double tcov = stod(fields[3]);

      if (tid <= 0 || qcov < matepar.overlap || tcov < matepar.overlap)
        continue;
      if ((size_t) tid > qry_seqs) {
        say("2sss$",
            "Warning: parsed query id out of range:",
            fields[0].c_str(), " - ignored");
        continue;
      }

      result[tid - 1]++;
    } catch (const exception &e) {
      // Parsing error - skip line
      continue;
    }
  }
  filein.close();

  remove(fasta);
  remove(db);
  remove(out);
  remove(log);
}

} // namespace

size_t nonpareil_mate_usearch(
      int *&result, char *file, char *sampleFile, int threads,
      size_t qry_seqs, unsigned int total_seqs, matepar_t matepar) {
  std::vector<shard_t> shards;
  size_t n_shards = 0;

  if (processID == 0) {
    int    co_located = ranks_on_this_node();
    double ram_mb_per_rank = matepar.ram_max_mb / (double) co_located;
    say(
      "4sisfs$", "RAM budget: ", co_located,
      " rank(s) sharing this node, ", ram_mb_per_rank, " Mb/rank"
    );

    shards = plan_shards(file, total_seqs, ram_mb_per_rank);
    n_shards = shards.size();
    say("3sus$", "Sharded subject dataset into ", (unsigned int) n_shards,
        " shard(s)");
  }

  // Broadcast the shard plan to all ranks
  broadcast_size_t(&n_shards);
  if (processID != 0) shards.resize(n_shards);
  for (size_t s = 0; s < n_shards; s++) {
    broadcast_size_t(&shards[s].start);
    broadcast_int(&shards[s].count);
    broadcast_int(&shards[s].max_len);
    broadcast_size_t(&shards[s].total_bases);
  }

  // Write the shards owned by this rank, then process them round-robin
  write_owned_shards_to_fasta(file, shards);
  for (size_t s = 0; s < n_shards; s++) {
    if ((int)(s % processes) != processID) continue;
    if (processID == 0)
      say("6susu$", "Processing shard ", (unsigned int)(s + 1), "/",
          (unsigned int) n_shards);
    process_one_shard(
      result, sampleFile, shards[s], threads, qry_seqs, matepar
    );
  }

  // Reduce multi-node results
  barrier_multinode();
  if (processes > 1) {
    int *result_sum = new int[qry_seqs];
    reduce_sum_int(result, result_sum, qry_seqs);
    if (processID == 0)
      for (size_t a = 0; a < qry_seqs; a++) result[a] = result_sum[a];
    delete[] result_sum;
  }
  barrier_multinode();

  return qry_seqs;
}

size_t nonpareil_mate(
    int *&result, char *file, int threads, unsigned int lines_in_ram,
    unsigned int total_seqs, unsigned int largest_seq, matepar_t matepar) {
  // Use `file` for both query and subject sequences
  return nonpareil_mate(
    result, file, file, threads, lines_in_ram, total_seqs, largest_seq,
    largest_seq, matepar
  );
}

size_t nonpareil_mate(
      int *&result, char *file, char *q_file, int threads,
      unsigned int lines_in_ram, unsigned int total_seqs,
      unsigned int largest_seq, unsigned int q_largest_seq, matepar_t matepar) {
  // Vars
  int    no_blocks_sbj = 0, no_blocks_qry = 0, no_seqs_block_qry = 0,
         no_seqs_block_sbj = 0, tmp_ram, result_i = 0, size_blockA, size_blockB;
  size_t qry_seqs = 0;
  char   **blockA, **blockB, *sampleFile;

  // Set subsampling
  sampleFile = new char[LARGEST_PATH];
  if (processID == 0) {
    std::filesystem::path tmp_path = tmp_dir();
    snprintf(sampleFile, LARGEST_PATH, "%s/subsample", tmp_path.c_str());
    say("3ss$", "Building query set at ", sampleFile);
    qry_seqs = sub_sample_seqs(
      q_file, sampleFile, matepar.qryportion, (char *)"enveomics-seq"
    );
    say("4sus$", "Query set built with ", qry_seqs, " sequences");
    if (qry_seqs == 0)
      error(
        "Impossible to create the query set.  Is the -X/-x value too small?");
  }
  broadcast_char(sampleFile, LARGEST_PATH);
  broadcast_int(&qry_seqs);

  // Blank results
  result = new int[qry_seqs];
  for (size_t a = 0; a < qry_seqs; a++) result[a] = 0;

  // If `-T usearch`, delegate to the specialized function
  if (matepar.type == 3) {
    return nonpareil_mate_usearch(
      result, file, sampleFile, threads, qry_seqs, total_seqs, matepar
    );
  }

  // Design blocks
  if (processID == 0){
    say("5sis$", "Designing the blocks scheme for ", total_seqs, " sequences");

    no_blocks_qry = (int)ceil(
      (double)qry_seqs * 2 / (double)lines_in_ram
    ); // Maximum half of the available slots
    if (no_blocks_qry == 0) no_blocks_qry = 1; // <-- Because of float precision
    no_seqs_block_qry = (int)ceil((double)qry_seqs / (double)no_blocks_qry);
    say("6sisi$",
        "Qry blocks:", no_blocks_qry, ", seqs/block:", no_seqs_block_qry);

    no_blocks_sbj = (int)ceil(
      (double)total_seqs / (double)(lines_in_ram - no_seqs_block_qry)
    );
    if (no_blocks_sbj == 0) no_blocks_sbj = 1; // <-- Because of float precision
    no_blocks_sbj = (int)ceil(
      (double)no_blocks_sbj / (double)processes
    ) * processes;
    no_seqs_block_sbj = (int)ceil((double)total_seqs / (double)no_blocks_sbj);
    say("6sisi$",
        "Sbj blocks:", no_blocks_sbj, ", seqs/block:", no_seqs_block_sbj);
  }
  broadcast_int(&no_blocks_qry);
  broadcast_int(&no_seqs_block_qry);
  broadcast_int(&no_blocks_sbj);
  broadcast_int(&no_seqs_block_sbj);

  // Mating
  if (processID == 0)
    say("3sisis$",
        "Mating sequences in ", no_blocks_qry,
        " by ", no_blocks_sbj, " blocks");
  if (processID == 0 && processes > 1)
    say("3sis$", "Silencing log in worker processes (", processes - 1, ")");
  for (int i = 0; i < no_blocks_qry; i++) {
    // Sequences in block A (qry)
    tmp_ram = (int) (
      ((double)no_seqs_block_qry / 1024) *
        q_largest_seq * (sizeof **blockA) / 1024
    );
    if (processID == 0)
      say("5sisi$",
          "Allocating ~", tmp_ram, " Mib in RAM for block qry:", i + 1);
    size_blockA = get_seqs(
      blockA, sampleFile, i * no_seqs_block_qry + 1, no_seqs_block_qry,
      q_largest_seq, (char *)"enveomics-seq"
    );
    if (size_blockA == 0) error("Impossible to get the i-th query block", i);

    // Sequences in block B (sbj)
    for (int j = 0; j < no_blocks_sbj; j++) {
      if (j % processes == processID) {
        tmp_ram = (int)(
          ((double)no_seqs_block_sbj / 1024) *
            largest_seq * (sizeof **blockB) / 1024
        );
        if (processID == 0)
          say("5sisi$",
              "Allocating ~", tmp_ram, " Mib in RAM for block sbj:", j + 1);
        size_blockB = get_seqs(
          blockB, file, j * no_seqs_block_sbj + 1, no_seqs_block_sbj,
          largest_seq, (char *)"enveomics-seq"
        );
        if (size_blockB == 0)
          error("Impossible to get the i-th subject block", j);

        // Mate
        if (processID == 0)
          say("4sisi$",
              "Computing block ", (i + 1) * (j + 1),
              "/", no_blocks_qry * no_blocks_sbj);
        nonpareil_count_mates_block(
          result, result_i, blockA, blockB, size_blockA, size_blockB, threads,
          matepar
        );
        for (int a = 0; a < size_blockB; a++) delete [] blockB[a];
        delete[] blockB;
      }
    }
    result_i += size_blockA;
    for (int a = 0; a < size_blockA; a++) delete [] blockA[a];
    delete[] blockA;
  }
  barrier_multinode();

  // Reduce multi-node results
  if (processes > 1) {
    int *result_sum = new int[qry_seqs];
    reduce_sum_int(result, result_sum, qry_seqs);
    if (processID == 0)
      for (size_t a = 0; a < qry_seqs; a++) result[a] = result_sum[a];
  }
  barrier_multinode();

  if (processID == 0) remove(sampleFile);

  return qry_seqs;
}

void nonpareil_count_mates_block(
      int *&result, int from_in_result, char **&blockA, char **&blockB,
      int sizeBlockA, int sizeBlockB, int threads, matepar_t matepar) {
  // Vars
  if (sizeBlockA < threads) threads = sizeBlockA;
  pthread_t       thread[threads];
  matejob_t       matejob[threads];
  pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
  unsigned int    mates_per_thr = (unsigned int) ceil(
                        (double) sizeBlockA / threads);
  int             rc;

  // Set threads
  threads = (int) ceil((double) sizeBlockA / mates_per_thr);

  // Launch jobs
  if (processID == 0)
    say("5sis$",
        "Launching parallel comparisons to ", threads, " threads");
  for (int thr = 0; thr < threads; thr++) {
    matejob[thr].id = thr; // The ID of the thread
    matejob[thr].from = mates_per_thr * thr; // The first qry sequence to
                                             // process (in zero-count)
    matejob[thr].number = (
      matejob[thr].from + mates_per_thr > (size_t) sizeBlockA ?
        sizeBlockA - matejob[thr].from : mates_per_thr
    ); // How many qry sequences to process
    matejob[thr].from_in_result = from_in_result; // Where to start saving
                                                  // results (in zero-count)
    matejob[thr].par = matepar; // It's cheap to create multiple copies of this,
                                // and it's safer than passing a reference.
    matejob[thr].result = &result; // This is only used with mutex (to save RAM)
    matejob[thr].mutex = &mutex; // And this is the mutex, it MUST be the same
    matejob[thr].blockA = &blockA;
    matejob[thr].blockB = &blockB;
    matejob[thr].size_blockA = sizeBlockA;
    matejob[thr].size_blockB = sizeBlockB;

    if (matejob[thr].number == 0)
       error("Unexpectedly, the thread contains zero load", thr);

    if ((rc = pthread_create(
          &thread[thr], NULL, &nonpareil_count_mates_thr, (void *) &matejob[thr]
        ))) error("Thread creation failed", (char) rc);
  }

  // Gather jobs
  for (int thr = 0; thr < threads; thr++) pthread_join(thread[thr], NULL);

  return;
}

void *nonpareil_count_mates_thr(void *matejob_ref) {
  // Vars
  matejob_t *matejob = (matejob_t *) matejob_ref;
  int       *result_cp = new int[matejob->number];

  if (!result_cp)
    error("Impossible to allocate memory for the results of the new thread",
          matejob->id);

  // Run comparisons
  for (size_t a = 0; a < matejob->number; a++) result_cp[a] = 0;
  nonpareil_count_mates(
    result_cp, *matejob->blockA, *matejob->blockB, matejob->from,
    matejob->number, 0, matejob->size_blockB,
    (matejob->id == 0 ? (int) ceil((double) matejob->number / 100.0) : 0),
    matejob->par
  );

  // Transfer the results to the external array
  pthread_mutex_lock( matejob->mutex );
    if (processID == 0)
      say("4sisis>", "Thread ", matejob->id, " completed ", matejob->number,
          " comparisons, joining results");
    int *&result_ref = *matejob->result;

    // result_ref[position + first of the block + first of the thread]
    for (size_t i = 0; i < matejob->number; i++)
      result_ref[i + matejob->from_in_result + matejob->from] += result_cp[i];
   pthread_mutex_unlock( matejob->mutex );

   return (void *) 0;
}

void nonpareil_count_mates(
      int *&result, char **&blockA, char **&blockB, int fromA, int numberA,
      int fromB, int numberB, int talk, matepar_t matepar) {
  // Finally, this is the core of the per-block-per-thread comparisons
  for (int i = 0; i < numberA; i++) {
    if (processID == 0 && talk > 0 && i % talk == 0)
      say("4sfs^", "Searching sequences: ", (double) i * 100 / (double) numberA,
          "% of the block");
    for (int j = 0; j < numberB; j++)
      if (nonpareil_compare_reads(
            blockA[i + fromA], blockB[j + fromB], matepar))
        result[i]++;
  }

  return;
}

bool nonpareil_compare_reads(char *seqA, char *seqB, matepar_t matepar) {
  // Vars
  int lenA = strlen(seqA), lenB = strlen(seqB);

  // Compare
  if (lenA < lenB)
    return nonpareil_compare_reads_shortfirst(seqA, seqB, lenA, lenB, matepar);
  return nonpareil_compare_reads_shortfirst(seqB, seqA, lenB, lenA, matepar);
}

bool nonpareil_compare_reads_shortfirst(
      char *seqA, char *seqB, int lenA, int lenB, matepar_t matepar) {
  int    min_len = (int) ceil(matepar.overlap * lenA);
  double disimilarity = 1.0 - matepar.similarity;

  // Compare strands W-W
  for (int i = min_len - lenB; i <= lenA-min_len; i++) {
    int errors = 0, from_seqA = (i < 0 ? 0 : i),
        len_within = (i < 0 ? lenB+i : lenA-i),
        max_errors = (int) (disimilarity * (double) len_within),
        from_seqB = (i < 0 ? -1 * i : 0);

    for (int j = 0; j < len_within; j++) {
       if (
             (seqA[j + from_seqA] != seqB[j + from_seqB] ||
               (matepar.n_as_mismatch &&
                 (seqA[j + from_seqA] == 'N' ||
                   seqB[j + from_seqB] == 'N'))) &&
             (++errors > max_errors))
          goto next_i; // Cannot be 'continue', because it's a nested loop
    }
    return true;

    next_i: ; // To avoid problems with the nested loop.
  }

  // Compare strands C-W
  if (matepar.revcom) {
    char *rc = new char[lenA + 1];
    bool out;

    reverse_complement(rc, seqA);
    matepar.revcom = false;
    out = nonpareil_compare_reads_shortfirst(rc, seqB, lenA, lenB, matepar);
    matepar.revcom = true; // Not really necesary, just in case it's passed as
                           // reference or something like that

    free(rc);
    return out;
  }

  return false;
}

void nonpareil_save_mates(int *&result, int no_results, char *file) {
   // Vars
   ofstream fileh;

   fileh.open(file, ios::out);
   if (!fileh.is_open()) error("Cannot open the file", file);
   for (int a = 0; a < no_results; a++) fileh << result[a] << endl;

   return;
}
