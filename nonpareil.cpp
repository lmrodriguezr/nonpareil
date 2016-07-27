// nonpareil - Calculation of nonpareil curves
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @version 2.4
// @license artistic license 2.0

/*
  ROADMAP

  o To validate "logarithmic sampling".

*/

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
#define NP_VERSION 3.000

using namespace std;
int processID;
int processes;

void help(const char *msg){
   if(processID==0 && msg!=NULL && strlen(msg) != 0) cerr << endl << msg << endl << endl;
   if(processID==0) cerr << "DESCRIPTION"<< endl
   	<< "   Nonpareil uses the redundancy of the reads in metagenomic datasets to estimate the average coverage and predict the" << endl
	<< "   amount of sequences that will be required to achieve 'nearly complete coverage'." << endl
	<< endl
   	<< "USAGE" << endl
  << "   nonpareil -s sequences.fa -T alignment -b output [options]" << endl
	<< "   nonpareil -s sequences.fa -T kmer -f fastq -b output [options]" << endl
	<< "   nonpareil -h" << endl
	<< "   nonpareil -V" << endl
	<< endl
   	<< "MANDATORY ARGUMENTS" << endl
	<< "   -s <str> : Path to the (input) file containing the sequences.  This is lowercase S." << endl
	<< endl
	<< "COMMON OPTIONS" << endl
	<< "   -f <str> : The format of the sequences.  Can be 'fasta' or 'fastq'.  By default: 'fasta'." << endl
	<< "   -b <str> : Path to the prefix for all the output files.  Replaces the options: -a, -C, -l, and -o; generating files" << endl
	<< "              with the suffixes .npa, npc, .npl, and .npo, respectively, unless explicitly set." << endl
	<< "   -d <num> : Subsample iteratively applying this factor to the number of reads, resulting in logarithmic subsampling." << endl
	<< "              Use -d 0 to fall back to linear sampling, controlled by -m, -M, & -i (this was the default before v2.4)." << endl
	<< "              By default: 0.7." << endl
	<< "   -n <int> : Number of sub-samples to generate per point.  If it is not a multiple of the number of threads (see -t)," << endl
	<< "              it is rounded to the next (upper) multiple.  By default: 1024." << endl
	<< "   -L <num> : Minimum overlapping percentage of the aligned region on the largest sequence. The similarity (see -S) is" << endl
	<< "              evaluated for the aligned region only.  By default: 50." << endl
	<< "   -X <int> : Maximum number of reads to use as query.  This is capital X.  By default, 1,000 reads." << endl
	<< "   -q <str> : Path to the (input) file containing a second dataset to be used as query, for dataset comparisons.  This" << endl
	<< "              option is currently experimental." << endl
	<< "   -R <int> : Maximum RAM usage in Mib.  Ideally this value should be larger than the sequences to analyze (discarding" << endl
	<< "              non-sequence elements like headers or quality).  This is particularly important when running in multiple" << endl
	<< "              cores (see -t).  This value is approximated.  By default 1024." << endl
	<< "              Maximum value in this version: " << (UINT_MAX/1024) << endl
	<< "   -t <int> : Number of threads.  Highest efficiency when the number of sub-samples (see -n) is multiple of the number" << endl
	<< "              of threads.  By default: 2." << endl
	<< "   -v <int> : Verbosity level, for debugging purposes.  By default 7.  This is lowercase V." << endl
	<< "   -V       : Show version information and exit.  This is uppercase V." << endl
	<< "   -h       : Display this message and exit." << endl
  << "   -k <int> : kmer length. By default: 24" << endl
  << "   -T <str> : Nonpareil algorithm kmer or alignment accepted"
	<< endl
	<< "See all supported arguments and additional documentation at http://nonpareil.readthedocs.org or execute man nonpareil." << endl
	<< endl;
   finalize_multinode();
   if(processID==0){ exit(1); }else{ exit(0); }
}

int main(int argc, char *argv[]) {
   init_multinode(argc, argv, processID, processes);
   if(processID==0) cout << "Nonpareil v" << NP_VERSION << endl;
   if(argc<=1) help("");

   // Vars
   char		*file, *format=(char *)"fasta", *nonpareiltype, *alldata, *cntfile, *outfile, *namFile,
   		*seqFile, *baseout, *qfile, *qNamFile, *qSeqFile;
   double	min=0.0, max=1.0, itv=0.01, qry_portion=0, min_sim=0.95, ovl=0.50, *sample_result,
   		avg_seq_len, adj_avg_seq_len, divide=0.7, q_avg_seq_len;
   int		v=7, largest_seq, rseed=time(NULL), n=1024, k=24, thr=2, ram=1024, *mates, samples_no,
   		sample_i, sample_after_20, sampling_points, q_largest_seq;
   unsigned int	q_total_seqs, lines_in_ram, hX=1000, qry_seqs_no, ram_Kb, required_ram_Kb;
   unsigned long long int total_seqs;
   bool		n_as_mismatch=false, portion_label=false, revcom=true, ok, autoadjust=false,
   		alt_query=false;
   matepar_t	matepar;
   samplepar_t 	samplepar;

   alldata = (char *)"";
   cntfile = (char *)"";
   outfile = (char *)"-";


   // GetOpt
   int		optchr;
   while ((optchr = getopt (argc, argv, "Aa:b:cC:d:f:Fhi:k:l:L:m:M:n:No:q:r:R:s:S:t:T:v:Vx:X:")) != EOF)
      switch(optchr) {
	 case 'a': alldata = optarg;		break;
	 case 'A': autoadjust = true;		break;
	 case 'b': baseout = optarg;		break;
	 case 'c': revcom = false;		break;
	 case 'C': cntfile = optarg;		break;
	 case 'd': divide = atof(optarg);	break;
	 case 'f': format = optarg;		break;
	 case 'F': portion_label = true;	break;
	 case 'h': help("");			break;
	 case 'i': itv = atof(optarg);		break;
   case 'k': k = atoi(optarg);  break;
	 case 'l': open_log(optarg);		break;
	 case 'L': ovl=atof(optarg)/100.0;	break;
	 case 'm': min = atof(optarg);		break;
	 case 'M': max = atof(optarg);		break;
	 case 'n': n = atoi(optarg);		break;
	 case 'N': n_as_mismatch=true;		break;
	 case 'o': outfile = optarg;		break;
	 case 'q': qfile = optarg; alt_query = true;	break;
	 case 'r': rseed=atoi(optarg);		break;
	 case 'R': ram = atoi(optarg);		break;
   case 's': file = optarg;		break;
	 case 'S': min_sim=atof(optarg);	break;
	 case 't': thr = atoi(optarg);		break;
	 case 'T': nonpareiltype = optarg; break; // For backwards compatibility
	 case 'v': v = atoi(optarg);		break;
	 case 'V': finalize_multinode(); return 0;
	 case 'x': qry_portion = atof(optarg);	break;
	 case 'X': hX = atoi(optarg);		break;
     }
   // Initialize
   set_verbosity(v);
   if(strlen(nonpareiltype)==0) help("");
   if(strcmp(nonpareiltype,"kmer")!=0 && strcmp(nonpareiltype,"alignment")!=0)
        help("Bad argument for -t option, accepted values are kmer or alignment");
   if(strlen(file)==0) help("");
   if(strcmp(format, "fasta")!=0 & strcmp(format, "fastq")!=0)
   				help("Unsupported value for -f option");
   if(min<0 | min>1)		help("Bad argument for -m option, accepted values are numbers in the range [0, 1]");
   if(max<0 | max>1)		help("Bad argument for -M option, accepted values are numbers in the range [0, 1]");
   if(itv<=0 | itv>1)		help("Bad argument for -i option, accepted values are numbers in the range (0, 1]");
   if(ovl<=0.0 | ovl>1.0)	help("Bad argument for -L option, accepted values are numbers in the range (0, 100]");
   if(thr<=0)			help("Bad argumement for -t option, accepted values are positive non-zero integers");
   if(n<=0)			help("Bad argument for -n option, accepted values are positive non-zero integers");
   if(ram<=0)			help("Bad argument for -R option, accepted values are positive non-zero integers");
   if(min_sim<=0 | min_sim>1)	help("Bad argument for -S option, accepted values are numbers in the range (0, 1]");
   if(thr<=0)			help("Bad argument for -t option, accepted values are positive non-zero integers");
   if(qry_portion<0 | qry_portion>1)
   				help("Bad argument for -x option, accepted values are numbers in the range (0, 1]");
   if(divide<0 | divide>=1)	help("Bad argument for -d option, accepted values are numbers in the range (0, 1)");
   char alldataTmp[LARGEST_PATH], outfileTmp[LARGEST_PATH], cntfileTmp[LARGEST_PATH];
   if(baseout && (strlen(baseout)>0)){
      if(!alldata || (strlen(alldata)<=0)){ sprintf(alldataTmp, "%s.npa", baseout); alldata=alldataTmp; }
      if(!cntfile || (strlen(cntfile)<=0)){ sprintf(cntfileTmp, "%s.npc", baseout); cntfile=cntfileTmp; }
      if(!outfile || (strcmp(outfile, "-")==0)){ sprintf(outfileTmp, "%s.npo", baseout); outfile=outfileTmp; }
      if(processID==0 && !log_is_open()) {
         char logfile[LARGEST_PATH];
	 sprintf(logfile, "%s.npl", baseout);
	 open_log(logfile);
      }
   }

   if(alldata && (strlen(alldata)>0)) remove(alldata);
   if(cntfile && (strlen(cntfile)>0)) remove(cntfile);
   if(outfile && (strlen(outfile)>0) & (strcmp(outfile, "-")!=0)) remove(outfile);

   if(strcmp(nonpareiltype,"kmer")==0) {
     if(alt_query) {
       if(strcmp(format,"fasta")==0) {
         cout << "testing" << endl;
         ifstream qifs((string(qfile)));
         say("1ss$","reading query", file);
         FastaReader qfastaReader(qifs);
         References references = References(qfastaReader, k, alt_query);
         say("1s$","Started counting");
         ifstream ifs((string(file)));
         FastaReader fastaReader(ifs);
         KmerCounter counter = KmerCounter(references, fastaReader, string(cntfile));
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
     else {
       if(strcmp(format,"fastq")==0) {
         ifstream ifs((string(file)));
         say("1ss$","reading ", file);
         FastqReader fastqReader(ifs);
         say("1sus$","Picking ", hX, " random sequences");
         References references = References(fastqReader, hX, k);
         say("1s$","Started counting");
         KmerCounter counter = KmerCounter(references, fastqReader, string(cntfile));
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
       else if(strcmp(format,"fasta")==0) {
         ifstream ifs((string(file)));
         say("1ss$","reading ", file);
         FastaReader fastaReader(ifs);
         say("1sus$","Picking ", hX, " random sequences");
         References references = References(fastaReader, hX, k);
         say("1s$","Started counting");
         KmerCounter counter = KmerCounter(references, fastaReader, string(cntfile));
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

   srand(rseed+processID);
   say("9si$", "Hello from worker ", processID);
   barrier_multinode();

   // Parse file
   if(processID==0) say("1s$", "Counting sequences");
   namFile = (char *)malloc(LARGEST_PATH * (sizeof *namFile));
   seqFile = (char *)malloc(LARGEST_PATH * (sizeof *seqFile));
   if(processID==0){
      total_seqs = build_index(file, format, namFile, seqFile, largest_seq, avg_seq_len);
      if(largest_seq<1) error("Your sequences are empty or an internal error occurred.  Largest sequence is ", largest_seq);
      say("2sss$", "The file ", seqFile, " was just created");
      say("4sis$", "Longest sequence has ", largest_seq, " characters");
      say("4sfs$", "Average read length is ", avg_seq_len, " bp");
      if(total_seqs==0) error("The file you provided does not contain sequences.  Before re-running please delete the file", seqFile);
      say("1sus$", "Reading file with ", total_seqs, " sequences");
      say("9s$", "Broadcasting");
   }
   total_seqs = broadcast_int(total_seqs);
   namFile = broadcast_char(namFile, LARGEST_PATH);
   seqFile = broadcast_char(seqFile, LARGEST_PATH);
   largest_seq = broadcast_int(largest_seq);
   avg_seq_len = broadcast_double(avg_seq_len);
   barrier_multinode();

   // Parse Q-File
   qNamFile = (char *)malloc(LARGEST_PATH * (sizeof *qNamFile));
   qSeqFile = (char *)malloc(LARGEST_PATH * (sizeof *qSeqFile));
   if(processID==0){
      say("1s$", "Counting query sequences");
      if(alt_query){
	 q_total_seqs = build_index(qfile, format, qNamFile, qSeqFile, q_largest_seq, q_avg_seq_len);
	 if(q_largest_seq<1) error("Your query sequences are empty or an internal error occurred.  Largest sequence is ", q_largest_seq);
	 say("2sss$", "The file ", qSeqFile, " was just created");
	 say("4sis$", "Longest query sequence has ", q_largest_seq, " characters");
	 say("4sfs$", "Average query read length is ", q_avg_seq_len, " bp");
	 if(q_total_seqs==0) error("The query file you provided does not contain sequences.  Before re-running please delete the file", qSeqFile);
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
   q_total_seqs = broadcast_int(q_total_seqs);
   qNamFile = broadcast_char(qNamFile, LARGEST_PATH);
   qSeqFile = broadcast_char(qSeqFile, LARGEST_PATH);
   q_largest_seq = broadcast_int(q_largest_seq);
   q_avg_seq_len = broadcast_double(q_avg_seq_len);
   barrier_multinode();


restart_vars:
   say("9sis$", "Worker ", processID, " @start_vars.");
   if(processID==0){
      // Re-wire query portion
      if(qry_portion!=0) hX = (size_t)total_seqs*qry_portion;
      qry_portion = (double)hX/(alt_query ? q_total_seqs : total_seqs);

      // Prepare memory arguments
      if(ram > UINT_MAX/1024) error("The memory to allocate is huge, I cannot manage such a big number.  Reduce -R and try again", ram);
      ram_Kb = ram*1024;
      required_ram_Kb = 2*(int)hX*sizeof(int)*thr/1024 + 2048; //<- 2Mb to store other things (usually <1Mb)
      if(ram_Kb < required_ram_Kb)
	 error("The amount of memory allowed is lesser than the minimum required, increase -R to at least", (double)required_ram_Kb/1024);
      lines_in_ram = (ram_Kb - required_ram_Kb)/(largest_seq + 3); // <- This is Kibi-lines, not lines
      if(lines_in_ram > UINT_MAX/1024){
	 say("1sfs$", "WARNING: Unable to represent RAM in bits, lowering to ", (double)UINT_MAX/(1024*1024), "Mb");
	 lines_in_ram = UINT_MAX;
      }else lines_in_ram *= 1024; // <- Beware, I am rounding down to the Kibi.  This is NOT the actual count.
      say("3sfsisfs$",
	   "Sequences to store in ", (double)(ram_Kb - required_ram_Kb)/1024, "Mb free: ", lines_in_ram,
	   " (", (double)lines_in_ram*100.0/total_seqs, "%)");
   }
   hX = broadcast_int(hX);
   qry_portion = broadcast_double(qry_portion);
   ram_Kb = broadcast_int(ram_Kb);
   required_ram_Kb = broadcast_int(required_ram_Kb);
   lines_in_ram = broadcast_int(lines_in_ram);
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
   if(processID==0) say("1sfsis$", "Querying library with ", qry_portion, " times the total size (", hX," seqs)");
   if(alt_query){
      qry_seqs_no = nonpareil_mate(mates, seqFile, qSeqFile, thr, lines_in_ram, total_seqs, largest_seq, q_largest_seq, matepar);
   }else{
      qry_seqs_no = nonpareil_mate(mates, seqFile, thr, lines_in_ram, total_seqs, largest_seq, matepar);
   }
   if(processID==0){
      if(alt_query) for(size_t a=0; a<qry_seqs_no; a++) mates[a]++; // <- Accounts for the lack of self-matches
      if(cntfile && (strlen(cntfile)>0)) nonpareil_save_mates(mates, qry_seqs_no, cntfile);
   }
   barrier_multinode();

   // Sampling
restart_samples:
   say("9sis$", "Worker ", processID, " @start_samples.");
   sampling_points = (divide==0) ? ((int)ceil((max-min)/itv)+1) : ((int)ceil( (log(2) - log(total_seqs))/log(divide) )+2);
   sample_t	sample_summary[sampling_points];
   size_t	dummy=0;
   if(processID==0) {
     if(strcmp(nonpareiltype,"kmer")!=0) {
      sample_i=sample_after_20=0;
      samplepar.np_version = NP_VERSION;
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
      samplepar.type = 1;
    }
    else {
      sample_i=sample_after_20=0;
      samplepar.np_version = NP_VERSION;
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
      samplepar.k = k;
      samplepar.adj_avg_read_len = adj_avg_seq_len;
      samplepar.type = 2; //kmer
    }

      say("1s$", "Sub-sampling library");
      double a=min;

      while(sample_i < sampling_points){
	 //if(sampling_points<=sample_i)
	   //error("Unexpected number of sampling points.  This can be due to a precision problem, try decreasing the -i parameter");
	 if(divide==0){
	    samplepar.portion = min + itv*sample_i;
	 }else{
	    samplepar.portion = sample_i==0 ? 0 : pow(divide, sampling_points-sample_i-1);
	 }
	 samplepar.replicates = n;

	 samples_no = nonpareil_sample_portion(sample_result, thr, samplepar);
	 //if(processID==0){
	    sample_summary[sample_i++] = nonpareil_sample_summary(sample_result, samples_no, alldata, outfile, samplepar);
	    if(samplepar.portion<=0.2) sample_after_20 = sample_i;
	 //}
      }
      dummy = 1;
   }
   dummy = broadcast_int(dummy);
   barrier_multinode();
   // Check results
restart_checkings:
   say("9sis$", "Worker ", processID, " @start_checkings.");
   if(processID==0){
      say("1s>", "Evaluating consistency");
      ok = true;
      // Low sequencing depth
      if(sample_after_20<=sample_i && sample_summary[sample_after_20].q2==0.0){
	 say("1s$", "WARNING: The estimation at 20% has median zero, possibly reflecting inaccurate estimations");
	 if(qry_portion<1.0 && hX<3000){
	    if(autoadjust){
	       hX = 0;
	       qry_portion *= 2.0;
	       if(qry_portion >= 1.0) qry_portion = 1.0;
	       say("1sf$", "AUTOADJUST: -x ", qry_portion);
	       goto restart_vars;
	    } else say("1sf$", "To increase the sensitivity of the estimations increase the -x parameter, currently set at ", qry_portion);
	 } else if(ovl > 0.25) {
	    if(autoadjust){
		    if(ovl>1.0)  ovl = 1.0;// This should never happen
	       else if(ovl>0.75) ovl = 0.75;
	       else if(ovl>0.5)  ovl = 0.5;
	       else if(ovl>0.25) ovl = 0.25;
	       else error("Impossible to reduce the -L parameter further, sequencing depth under detection level");
	       say("1sf$", "AUTOADJUST: -L ", ovl*100);
	       goto restart_mates;
	    }else say("1sf$", "To increase the sensitivity of the estimations decrease the -L parameter, currently set at ", ovl*100);
	 } else {
	    say("1s$", "The portion used as query (-x) is currently set to the maximum, and the overlap (-L) is set to the minimum");
	    say("1s$", "The dataset is probably too small for reliable estimations");
	    if(min_sim>0.75) say("1s$", "You could decrease the similarity (-S) but values other than 0.95 are untested");
	    error("Sequencing depth under detection limit.");
	 }
	 ok = false;
      }
      // High sequencing depth
      if(sample_summary[sample_i-1].avg >= 0.95){
	 say("1s$", "WARNING: The curve reached near-saturation, hence coverage estimations could be unreliable");
	 if(ovl<1.0){
	    if(autoadjust){
		    if(ovl<0.25)  ovl = 0.25;
	       else if(ovl<0.5) ovl = 0.5;
	       else if(ovl<0.75)  ovl = 0.75;
	       else if(ovl<1.0) ovl = 1.0;
	       say("1sf$", "AUTOADJUST: -L ", ovl*100.0);
	       goto restart_mates;
	    }else say("1sf$", "To avoid saturation increase the -L parameter, currently set at ", ovl*100);
	 } else {
	    say("1s$", "The overlap (-L) is currently set to the maximum, meaning that the actual coverage is probably above 100X");
	    if(min_sim<1.0) say("1s$", "You could increase the similarity (-S) but values other than 0.95 are untested");
	    error("Sequencing depth above detection limit.");
	 }
	 ok = false;
      }
      // Low resolution
      if(sample_i>5 && sample_summary[5].avg >= 0.95){
	 say("1s$", "WARNING: The curve reached near-saturation in 6 or less points, hence diversity estimations could be unreliable");
	 if(autoadjust){
	    itv *= 0.5;
	    say("1sf$", "AUTOADJUST: -i ", itv);
	    goto restart_samples;
	 }else say("1sf$", "To increase the resolution of the curve increase the -i parameter, currently set at ", itv);
	 ok = false;
      }
      if(ok) say("1s$", "Everything seems correct");
      dummy = 2;
   }
   dummy = broadcast_int(dummy);
   barrier_multinode();

exit:
   say("9sis$", "Worker ", processID, " @exit.");
   // Clean temporals
   if(processID==0){
      remove(namFile);
      remove(seqFile);
      close_log();
   }
   barrier_multinode();
   finalize_multinode();
   return 0;
}
