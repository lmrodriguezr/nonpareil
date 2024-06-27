// nonpareil - Calculation of nonpareil curves
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @author Santosh Kumar G.
// @license Artistic-2.0

#include <math.h>
#include <getopt.h>
#include "enveomics/universal.h"
#include "enveomics/multinode.h"
#include "enveomics/sequence.h"
#include "enveomics/nonpareil_mating.h"
#include "enveomics/nonpareil_sampling.h"
#include "enveomics/References.h"
#include "enveomics/KmerCounter.h"
#include <string>

#define LARGEST_PATH 4096
#define NP_VERSION "3.5.2"

using namespace std;
int processID;
int processes;

void help(const char *msg){
  if(processID == 0 && msg != NULL && strlen(msg) != 0)
    cerr << endl << msg << endl << endl;
  if(processID == 0)
    cerr
    <<"DESCRIPTION"                                                      << endl
    <<"  Nonpareil uses the redundancy of the reads in metagenomic"      << endl
    <<"  datasets to estimate the average coverage and predict the"      << endl
    <<"  amount of sequences that will be required to achieve 'nearly"   << endl
    <<"  complete coverage'."                                            << endl
    << endl
    <<"USAGE"                                                            << endl
    <<"  nonpareil -s sequences.fa -T alignment -b output [options]"     << endl
    <<"  nonpareil -s sequences.fq -T kmer -f fastq -b output [options]" << endl
    <<"  nonpareil -s seqs.fq.gz -T kmer -f fastq -b output [options]"   << endl
    <<"  nonpareil -h"                                                   << endl
    <<"  nonpareil -V"                                                   << endl
    << endl
    <<"MANDATORY ARGUMENTS"                                              << endl
    <<"  -s <str> : Path to the (input) file containing the sequences."  << endl
    <<"             Gzipped files are supported with .gz extension"      << endl
    <<"  -T <str> : Nonpareil algorithm, 'kmer' or 'alignment' accepted" << endl
    <<"  -f <str> : The format of the sequence: 'fasta' or 'fastq'"      << endl
    << endl
    <<"COMMON OPTIONS"                                                   << endl
    <<"  -b <str> : Path to the prefix for all the output files"         << endl
    <<"  -X <int> : Maximum number of reads to use as query"             << endl
    <<"             By default: 1000 for alignment, 10000 for kmer"      << endl
    <<"  -k <int> : kmer length. By default: 24"                         << endl
    <<"  -n <int> : Number of sub-samples to generate per point."        << endl
    <<"             If it is not a multiple of the number of threads"    << endl
    <<"             (-t), it is rounded to the next (upper) multiple."   << endl
    <<"             By default: 1024"                                    << endl
    <<"  -L <num> : Minimum overlapping percentage of the aligned region"<< endl
    <<"             on the largest sequence. The similarity (see -S) is" << endl
    <<"             evaluated for the aligned region only."              << endl
    <<"             By default: 50"                                      << endl
    <<"  -R <int> : Maximum RAM usage in Mib. Ideally this value should" << endl
    <<"             be larger than the sequences to analyze (discarding" << endl
    <<"             non-sequence elements like headers or quality). This"<< endl
    <<"             value is approximated. By default 1024"              << endl
    <<"             Maximum value in this version: " << (UINT_MAX/1024)  << endl
    <<"  -t <int> : Number of threads. By default: 2"                    << endl
    <<"  -v <int> : Verbosity level. By default 7"                       << endl
    <<"  -r <int> : Random seed to make runs reproducible. Currently"    << endl
    <<"             only implemented when -T alignment"                  << endl
    <<"  -V       : Show version information and exit"                   << endl
    <<"  -h       : Display this message and exit"                       << endl
    << endl
    <<"See all supported arguments and additional documentation at"      << endl
    <<"http://nonpareil.readthedocs.org or with `man nonpareil`"         << endl
    << endl;
  finalize_multinode();
  if (processID == 0) { exit(1); } else { exit(0); }
}

int main(int argc, char *argv[]) {
  init_multinode(argc, argv, processID, processes);
  if (processID == 0) cout << "Nonpareil v" << NP_VERSION << endl;
  if (argc <= 1) help("");

  // Vars
  char  *file = new char[LARGEST_PATH], *inputfile,
        *format=(char *)"", *nonpareiltype=(char *)"",
        *alldata, *cntfile, *outfile, *namFile, *seqFile, *baseout,
        *qfile, *qNamFile, *qSeqFile;
  double
        min=0.0, max=1.0, itv=0.01, qry_portion=0, min_sim=0.95, ovl=0.50,
        *sample_result, avg_seq_len=0.0, adj_avg_seq_len, divide=0.7,
        q_avg_seq_len=0.0;
  int   v=7, largest_seq=0, rseed=time(NULL), n=1024, k=24, thr=2, ram=1024,
        *mates, samples_no, sample_i, sample_after_20, sampling_points,
        q_largest_seq=0;
  unsigned int
        q_total_seqs=0, lines_in_ram, hX=0, qry_seqs_no, ram_Kb,
        required_ram_Kb;
  unsigned long long int total_seqs=0;
  bool  n_as_mismatch=false, portion_label=false, revcom=true, ok,
        autoadjust=false, alt_query=false, remove_input=false;
  matepar_t matepar;
  samplepar_t samplepar;

  alldata = (char *)"";
  cntfile = (char *)"";
  outfile = (char *)"-";


  // GetOpt
  int   optchr;
  while ((optchr = getopt (argc, argv,
          "Aa:b:cC:d:f:Fhi:k:l:L:m:M:n:No:q:r:R:s:S:t:T:v:Vx:X:")) != EOF) {
    switch (optchr) {
      case 'a': alldata = optarg;       break;
      case 'A': autoadjust = true;      break;
      case 'b': baseout = optarg;       break;
      case 'c': revcom = false;         break;
      case 'C': cntfile = optarg;       break;
      case 'd': divide = atof(optarg);  break;
      case 'f': format = optarg;        break;
      case 'F': portion_label = true;   break;
      case 'h': help("");               break;
      case 'i': itv = atof(optarg);     break;
      case 'k': k = atoi(optarg);       break;
      case 'l': open_log(optarg);       break;
      case 'L': ovl=atof(optarg)/100.0; break;
      case 'm': min = atof(optarg);     break;
      case 'M': max = atof(optarg);     break;
      case 'n': n = atoi(optarg);       break;
      case 'N': n_as_mismatch=true;     break;
      case 'o': outfile = optarg;       break;
      case 'q':
        qfile = optarg;
        alt_query = true;               break;
      case 'r': rseed=atoi(optarg);     break;
      case 'R': ram = (int)atoi(optarg);break;
      case 's': inputfile = optarg;     break;
      case 'S': min_sim=atof(optarg);   break;
      case 't': thr = atoi(optarg);     break;
      case 'T': nonpareiltype = optarg; break;
      case 'v': v = atoi(optarg);       break;
      case 'V':
        finalize_multinode(); return 0;
      case 'x':
        qry_portion = atof(optarg);     break;
      case 'X': hX = atoi(optarg);      break;
      default:
        help("Unrecognized flag");      break;
    }
  }
  // Set number of reads to use as query
  if (hX != 0) {
    // User-provided, do nothing
  } else if (strcmp(nonpareiltype,"kmer") == 0) {
    hX = 10000;
  } else if(strcmp(nonpareiltype,"alignment") == 0) {
    hX = 1000;
  }

  set_verbosity(v);
  if (strlen(nonpareiltype) == 0) help("-T is mandatory");
  if (strcmp(nonpareiltype, "kmer") != 0 &&
      strcmp(nonpareiltype, "alignment") != 0)
    help("Bad argument for -T option, accepted values are kmer or alignment");
  if (strlen(inputfile) == 0) help("-s is mandatory");
  if (strlen(format) == 0) help("-f is mandatory");
  if ((strcmp(format, "fasta") != 0) && (strcmp(format, "fastq") != 0))
    help("Unsupported value for -f option");
  if ((min < 0) | (min > 1))
    help("Bad argument for -m option, accepted range: [0, 1]");
  if ((max < 0) | (max > 1))
    help("Bad argument for -M option, accepted range: [0, 1]");
  if ((itv <= 0) | (itv > 1))
    help("Bad argument for -i option, accepted range: (0, 1]");
  if ((ovl <= 0.0) | (ovl > 1.0))
    help("Bad argument for -L option, accepted range: (0, 100]");
  if (thr <= 0)
    help("Bad argumement for -t option, accepted: positive non-zero integers");
  if (n <= 0)
    help("Bad argument for -n option, accepted: positive non-zero integers");
  if (ram <= 0)
    help("Bad argument for -R option, accepted: positive non-zero integers");
  if (thr <= 0)
    help("Bad argument for -t option, accepted: positive non-zero integers");
  if ((min_sim <= 0) | (min_sim > 1))
    help("Bad argument for -S option, accepted range: (0, 1]");
  if ((qry_portion < 0) | (qry_portion > 1))
    help("Bad argument for -x option, accepted range: (0, 1]");
  if ((divide < 0) | (divide >= 1))
    help("Bad argument for -d option, accepted range: (0, 1)");
  if ((k < 1) | (k > 32))
    help("Bad argument for -k option, accepted range: [1, 32]");

  char *alldataTmp = new char[LARGEST_PATH],
       *outfileTmp = new char[LARGEST_PATH],
       *cntfileTmp = new char[LARGEST_PATH];
  if (baseout && (strlen(baseout) > 0)) {
    if (!alldata || (strlen(alldata) <= 0)) {
      snprintf(alldataTmp, LARGEST_PATH, "%s.npa", baseout);
      alldata = alldataTmp;
    }
    if (!cntfile || (strlen(cntfile) <= 0)) {
      snprintf(cntfileTmp, LARGEST_PATH, "%s.npc", baseout);
      cntfile = cntfileTmp;
    }
    if (!outfile || (strcmp(outfile, "-") == 0)){
      snprintf(outfileTmp, LARGEST_PATH, "%s.npo", baseout);
      outfile = outfileTmp;
    }
    if (processID == 0 && !log_is_open()) {
      char logfile[LARGEST_PATH];
      snprintf(logfile, LARGEST_PATH, "%s.npl", baseout);
      open_log(logfile);
    }
  }
  broadcast_int(&rseed);
  barrier_multinode();
  srand(rseed + processID);

  // file checking
  int count = 0;
  int limit = 10 * hX; //metagenome should have 10 times more than query reads
  if (strcmp(nonpareiltype, "alignment") == 0) limit = hX;

  if (has_gz_ext(inputfile)) {
    remove_input = true;
    if (processID == 0) {
      snprintf(file, LARGEST_PATH, "%s.enve-tmp.%d", inputfile, getpid());
      gunz_file(inputfile, file);
      say("2sss$", "The file ", file, " was created");
    }
  } else {
    snprintf(file, LARGEST_PATH, "%s", inputfile);
  }
  broadcast_char(file, LARGEST_PATH);

  if (processID == 0) {
    Sequence test_temp;
    ifstream testifs((string(file)));
    if (strcmp(format, "fasta") == 0) {
      FastaReader testfastaReader(testifs);
      while(testfastaReader.readNextSeq(test_temp) != (size_t)(-1)) {
        count++;
        if (count > limit) break;
      }
    } else if(strcmp(format, "fastq") == 0) {
      FastqReader testfastqReader(testifs);
      while(testfastqReader.readNextSeq(test_temp) != (size_t)(-1)) {
        count++;
        if (count > limit) break;
      }
    } else {
      error("Unsupported format", format);
    }
  }
  broadcast_int(&count);
  if (count == 0) {
    error("No reads found, check that the input file exists and is readable");
  } else if (count < limit) {
    say("3sisis$", "Sequence count (", count, ") < limit (", limit, ")");
    if (strcmp(nonpareiltype, "alignment") == 0)
      error("Reduce the query reads (-X) to fit total reads", count);
    else
      error("Reduce the query reads (-X) to ≤ 10\% of total reads", count);
  }
  if(alldata && (strlen(alldata) > 0)) remove(alldata);
  if(cntfile && (strlen(cntfile) > 0)) remove(cntfile);
  if(outfile && (strlen(outfile) > 0) && (strcmp(outfile, "-") != 0))
    remove(outfile);
 
  if (strcmp(nonpareiltype, "kmer") == 0) {
    if (processID != 0) goto restart_samples;
    if (alt_query) {
      if(strcmp(format, "fasta") == 0) {
        say("1ss$", "WARNING: The kmer kernel implements an error correction ",
          "function only compatible with FastQ");
        ifstream qifs((string(qfile)));
        say("1ss$", "reading query", file);
        FastaReader qfastaReader(qifs);
        References references = References(qfastaReader, k, alt_query);
        say("1s$", "Started counting");
        ifstream ifs((string(file)));
        FastaReader fastaReader(ifs);
        KmerCounter counter = KmerCounter(references, fastaReader,
          string(cntfile));
        mates = new int[counter.getTotalSeqs()];
        counter.getCounts(mates);
        avg_seq_len = counter.getAvgLen();
        total_seqs = counter.getTotalSeqs();
        qry_seqs_no = counter.getTotalQSeqs();
        adj_avg_seq_len = avg_seq_len - k + 1;
        say("1sus$", "Read file with ", total_seqs, " sequences");
        say("1sfs$", "Average read length is ", avg_seq_len, "bp");
        goto restart_samples;
      }
    } else {
      if(strcmp(format, "fastq") == 0) {
        ifstream ifs((string(file)));
        say("1ss$", "reading ", file);
        FastqReader fastqReader(ifs);
        say("1sus$", "Picking ", hX, " random sequences");
        References references = References(fastqReader, hX, k);
        say("1s$", "Started counting");
        KmerCounter counter = KmerCounter(references, fastqReader,
          string(cntfile));
        mates = new int[hX];
        counter.getCounts(mates);
        avg_seq_len = counter.getAvgLen();
        total_seqs = counter.getTotalSeqs();
        qry_seqs_no = counter.getTotalQSeqs();
        adj_avg_seq_len = avg_seq_len - k + 1;
        say("1sus$", "Read file with ", total_seqs, " sequences");
        say("1sfs$", "Average read length is ", avg_seq_len, "bp");
        goto restart_samples;
      } else if(strcmp(format, "fasta") == 0) {
        say("1ss$","WARNING: The kmer kernel implements an error correction ",
          "function only compatible with FastQ");
        ifstream ifs((string(file)));
        say("1ss$", "reading ", file);
        FastaReader fastaReader(ifs);
        say("1sus$", "Picking ", hX, " random sequences");
        References references = References(fastaReader, hX, k);
        say("1s$", "Started counting");
        KmerCounter counter = KmerCounter(references, fastaReader,
          string(cntfile));
        mates = new int[hX];
        counter.getCounts(mates);
        avg_seq_len = counter.getAvgLen();
        total_seqs = counter.getTotalSeqs();
        qry_seqs_no = counter.getTotalQSeqs();
        adj_avg_seq_len = avg_seq_len - k + 1;
        say("1sus$", "Read file with ", total_seqs, " sequences");
        say("1sfs$", "Average read length is ", avg_seq_len, "bp");
        goto restart_samples;
      }
    }
  }

  say("9si$", "Hello from worker ", processID);
  barrier_multinode();

  // Parse file
  if(processID == 0) say("1s$", "Counting sequences");
  namFile = (char *)malloc(LARGEST_PATH * (sizeof *namFile));
  seqFile = (char *)malloc(LARGEST_PATH * (sizeof *seqFile));
  if(processID == 0){
    total_seqs = build_index(file, format, namFile, seqFile, largest_seq,
      avg_seq_len);
    if(largest_seq < 1)
      error("Empty sequences or internal error. Largest sequence: ",
        largest_seq);
    say("2sss$", "The file ", seqFile, " was just created");
    say("4sis$", "Longest sequence has ", largest_seq, " characters");
    say("4sfs$", "Average read length is ", avg_seq_len, " bp");
    if(total_seqs==0)
      error("No input sequences.  Before re-running please delete the file ",
        seqFile);
    say("1sus$", "Reading file with ", total_seqs, " sequences");
    say("9s$", "Broadcasting");
  }
  broadcast_int(&total_seqs);
  broadcast_char(namFile, LARGEST_PATH);
  broadcast_char(seqFile, LARGEST_PATH);
  broadcast_int(&largest_seq);
  broadcast_double(&avg_seq_len);
  barrier_multinode();

  // Parse Q-File
  qNamFile = (char *)malloc(LARGEST_PATH * (sizeof *qNamFile));
  qSeqFile = (char *)malloc(LARGEST_PATH * (sizeof *qSeqFile));
  if (processID == 0) {
    say("1s$", "Counting query sequences");
    if(alt_query){
      q_total_seqs = build_index(qfile, format, qNamFile, qSeqFile,
        q_largest_seq, q_avg_seq_len);
      if(q_largest_seq<1)
        error("No input sequences or internal error.  Largest sequence: ",
          q_largest_seq);
      say("2sss$", "The file ", qSeqFile, " was just created");
      say("4sis$", "Longest query sequence has ", q_largest_seq, " characters");
      say("4sfs$", "Average query read length is ", q_avg_seq_len, " bp");
      if(q_total_seqs==0)
        error("No input sequences.  Before re-running please delete the file ",
          qSeqFile);
      say("1sus$", "Reading query file with ", q_total_seqs, " sequences");
    }else{
      q_total_seqs=0;
      qNamFile=(char *)"";
      qSeqFile=(char *)"";
      q_largest_seq = 0;
      q_avg_seq_len = 0.0;
    }
    say("9s$", "Broadcasting");
  }
  broadcast_int(&q_total_seqs);
  broadcast_char(qNamFile, LARGEST_PATH);
  broadcast_char(qSeqFile, LARGEST_PATH);
  broadcast_int(&q_largest_seq);
  broadcast_double(&q_avg_seq_len);
  barrier_multinode();


restart_vars:
  say("9sis$", "Worker ", processID, " @start_vars.");
  if (processID == 0) {
    // Re-wire query portion
    if (qry_portion != 0) hX = (size_t)total_seqs*qry_portion;
    qry_portion = (double)hX/(alt_query ? q_total_seqs : total_seqs);

    // Prepare memory arguments
    if((size_t)ram > UINT_MAX/1024)
      error("The memory to allocate is too large, reduce -R", ram);
    ram_Kb = ram * 1024;
    required_ram_Kb = 2 * (int)hX * sizeof(int) * thr / 1024 + 2048;
    if(ram_Kb < required_ram_Kb)
      error("The amount of memory allowed is too small, increase -R to over ",
        (double)required_ram_Kb/1024);
    lines_in_ram = (ram_Kb - required_ram_Kb)/(largest_seq + 3);
    if(lines_in_ram > UINT_MAX/1024){
      say("1sfs$", "WARNING: Unable to represent RAM in bits, lowering to ",
        (double)UINT_MAX/(1024*1024), "Mb");
      lines_in_ram = UINT_MAX;
    }else lines_in_ram *= 1024; // <- Rounding down to the Kibi.
      say("3sfsisfs$",
        "Sequences to store in ", (double)(ram_Kb - required_ram_Kb)/1024,
        "Mb free: ", lines_in_ram, " (", (double)lines_in_ram*100.0/total_seqs,
        "%)");
  }
  broadcast_int(&hX);
  broadcast_double(&qry_portion);
  broadcast_int(&ram_Kb);
  broadcast_int(&required_ram_Kb);
  broadcast_int(&lines_in_ram);
  barrier_multinode();

  // Run comparisons
restart_mates:
  say("9sis$", "Worker ", processID, " @start_mates.");
  matepar.overlap = ovl;
  matepar.similarity = min_sim;
  matepar.qryportion = qry_portion;
  matepar.revcom = revcom;
  matepar.n_as_mismatch = n_as_mismatch;
  matepar.k = k;
  if (processID == 0)
    say("1sfsis$", "Querying library with ", qry_portion,
        " times the total size (", hX," seqs)");
  if (alt_query) {
    qry_seqs_no = nonpareil_mate(mates, seqFile, qSeqFile, thr, lines_in_ram,
      total_seqs, largest_seq, q_largest_seq, matepar);
  } else {
    qry_seqs_no = nonpareil_mate(mates, seqFile, thr, lines_in_ram, total_seqs,
      largest_seq, matepar);
  }
  if (processID == 0) {
    if (alt_query) for (size_t a=0; a<qry_seqs_no; a++) mates[a]++;
      // <- Accounts for the lack of self-matches
    if (cntfile && (strlen(cntfile) > 0))
      nonpareil_save_mates(mates, qry_seqs_no, cntfile);
  }
  barrier_multinode();

  // Sampling
restart_samples:
  say("9sis$", "Worker ", processID, " @start_samples.");
  sampling_points = (divide==0) ? ((int)ceil((max-min)/itv)+1) :
    ((int)ceil( (log(2) - log(total_seqs))/log(divide) )+2);
  sample_t sample_summary[sampling_points];
  size_t dummy=0;
  if (processID == 0) {
    sample_i = sample_after_20 = 0;
    samplepar.np_version = (char *)NP_VERSION;
    samplepar.replicates = n;
    samplepar.mates = &mates;
    samplepar.mates_size = qry_seqs_no;
    samplepar.portion_min = min;
    samplepar.portion_max = max;
    samplepar.portion_itv = itv;
    samplepar.seq_overlap = ovl;
    samplepar.total_reads = total_seqs;
    samplepar.max_read_len = largest_seq;
    samplepar.avg_read_len = avg_seq_len;
    samplepar.portion_as_label = portion_label;
    samplepar.divide = divide;

    if (strcmp(nonpareiltype, "kmer") != 0) {
      samplepar.type = 1;
    } else {
      samplepar.k = k;
      samplepar.adj_avg_read_len = adj_avg_seq_len;
      samplepar.type = 2; //kmer
    }

    say("1s$", "Sub-sampling library");
    while(sample_i < sampling_points){
      if (divide == 0) {
        samplepar.portion = min + itv*sample_i;
      } else {
        samplepar.portion = sample_i==0 ? 0 :
          pow(divide, sampling_points-sample_i-1);
      }
      samplepar.replicates = n;

      samples_no = nonpareil_sample_portion(sample_result, thr, samplepar);
      sample_summary[sample_i++] = nonpareil_sample_summary(sample_result,
        samples_no, alldata, outfile, samplepar);
      if (samplepar.portion <= 0.2) sample_after_20 = sample_i;
    }
    dummy = 1;
  }
  broadcast_int(&dummy);
  barrier_multinode();
  goto restart_checkings;

  // Check results
restart_checkings:
  say("9sis$", "Worker ", processID, " @start_checkings.");
  if(processID==0){
    say("1s>", "Evaluating consistency");
    ok = true;
    // Low sequencing depth
    if(sample_after_20<=sample_i && sample_summary[sample_after_20].q2==0.0){
      say("1ss$",
        "WARNING: The estimation at 20% has median zero, ",
        "possibly reflecting inaccurate estimations");
      if(qry_portion<1.0 && hX<3000){
        if(autoadjust){
          hX = 0;
          qry_portion *= 2.0;
          if(qry_portion >= 1.0) qry_portion = 1.0;
          say("1sf$", "AUTOADJUST: -x ", qry_portion);
          goto restart_vars;
        } else say("1sf$",
            "To increase the sensitivity increase -X, currently set at ", hX);
      } else if(ovl > 0.25) {
        if(autoadjust){
               if(ovl>1.0)  ovl = 1.0;// This should never happen
          else if(ovl>0.75) ovl = 0.75;
          else if(ovl>0.5)  ovl = 0.5;
          else if(ovl>0.25) ovl = 0.25;
          else error("Impossible to reduce -L further, sequencing depth under detection level");
          say("1sf$", "AUTOADJUST: -L ", ovl*100);
          goto restart_mates;
        } else say("1sf$",
            "To increase sensitivity, decrease -L, currently set at ", ovl*100);
      } else {
        say("1ss$",
          "The portion used as query (-x) is currently set to the maximum, ",
          "and the overlap (-L) is set to the minimum");
        say("1s$",
          "The dataset is probably too small for reliable estimations");
        if(min_sim>0.75) say("1s$",
          "You could decrease the -S but values other than 0.95 are untested");
        error("Sequencing depth under detection limit.");
      }
      ok = false;
    }
    // High sequencing depth
    if(sample_summary[sample_i-1].avg >= 0.95){
      if(ovl<1.0){
        if(autoadjust){
               if(ovl<0.25)  ovl = 0.25;
          else if(ovl<0.5) ovl = 0.5;
          else if(ovl<0.75)  ovl = 0.75;
          else if(ovl<1.0) ovl = 1.0;
          say("1sf$", "AUTOADJUST: -L ", ovl*100.0);
          goto restart_mates;
        }
      } else {
	  say("1ss$",
            "The overlap (-L) is currently set to the maximum, ",
            "meaning that the actual coverage is probably above 100X");
	  if(min_sim<1.0) say("1s$",
            "You could increase -S but values other than 0.95 are untested");
	  error("Sequencing depth above detection limit.");
      }
      ok = false;
    }
    // Low resolution
    if(sample_i>5 && sample_summary[5].avg >= 0.95){
      say("1ss$",
        "WARNING: The curve reached near-saturation in 6 or less points, ",
        "hence diversity estimations could be unreliable");
      if(autoadjust){
	itv *= 0.5;
	say("1sf$", "AUTOADJUST: -i ", itv);
	goto restart_samples;
      } else say("1ssf$",
        "To increase the resolution of the curve increase the -i parameter, ",
        "currently set at ", itv);
      ok = false;
    }
    if(ok) say("1s$", "Everything seems correct");
    dummy = 2;
  }
  broadcast_int(&dummy);
  barrier_multinode();
  goto exit;

exit:
  say("9sis$", "Worker ", processID, " @exit.");
  // Clean temporals
  if (processID == 0) {
    if (remove_input) remove(file);
    remove(namFile);
    remove(seqFile);
    close_log();
  }
  barrier_multinode();
  finalize_multinode();
  return 0;
}
