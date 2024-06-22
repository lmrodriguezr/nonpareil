// nuc_sampler - Estimates the composition of (poly-)nucleotides in a sequence sample
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @license artistic 2.0
// @version 2.0

#define _MULTI_THREADED
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <vector>
#include "enveomics/universal.h"
#include "enveomics/sequence.h"

//#define DEBUG(a) (cerr << "(LINE " << a << ")" << endl)
#define DEBUG(a) (a)
#define LARGEST_LINE	2048
#define LARGEST_PATH	2048
#define OMAG_STEP	1024

using namespace std;

typedef struct {
   char		*kmer;	// <- The polynucleotide
   unsigned int	label;	// <- A numeric representation of the kmer
   // Total count is calculated as: base*omag + hang
   unsigned int	base;	// <- The base in the above expression
   unsigned int	omag;	// <- The omag in the above expression
   unsigned int	hang;	// <- The hang in the above expression
} count_t;

int  V=0;
char nuc_char[] = {
  'N', 'A', 'C', 'G', 'T', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', '.'
};
char hex_char[] = {
  '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'
};

void help(const char *msg){
  if(msg!=NULL && strlen(msg) != 0) cerr << endl << msg << endl << endl;
  cerr
    << "DESCRIPTION"                                                     << endl
    <<"   Calculates or approximates the composition of nucleotides or"  << endl
    <<"   polynucleotides in large datasets"                             << endl
    << endl
    <<"USAGE"                                                            << endl
    <<"   nuc_sampler -s sequences.fa [options]"                         << endl
    <<"   nuc_sampler -h"                                                << endl
    <<"   nuc_sampler -V"                                                << endl
    << endl
    <<"MANDATORY ARGUMENTS"                                              << endl
    <<"   -s <str> : Path to the (input) file containing the sequences"  << endl
    << endl
    <<"ADDITIONAL OPTIONS"                                               << endl
    <<"   -o <str> : Path to the output file where results will be"      << endl
    <<"              saved. By default the results are sent to stdout."  << endl
    <<"              This is the same behavior as using a dash (-)."     << endl
    <<"              If an empty string is provided, does not produce"   << endl
    <<"              any output"                                         << endl
    <<"   -e       : Produce extended output"                            << endl
    <<"   -f <str> : The format of the sequences. Can be 'fasta' or"     << endl
    <<"              'fastq'. By default: 'fasta'"                       << endl
    <<"   -k <int> : Size of the polynucleotides (words) to count."      << endl
    <<"              By default: 1" << endl
    <<"   -r <int> : Random generator seed. By default current time"     << endl
    <<"   -v <int> : Verbose, for debugging purposes. By default 0"      << endl
    <<"   -x <num> : Probability of taking a sequence into account,"     << endl
    <<"              regardless of the sequence length. Lower values"    << endl
    <<"              reduce accuracy but increase speed. Any value"      << endl
    <<"              lower than 1 produces and approximation of the"     << endl
    <<"              composition. By default 1"                          << endl
    <<"   -h       : Display this message and exit"                      << endl
    <<"   -V       : Show version information and exit"                  << endl
    << endl
    <<"INPUT"                                                            << endl
    <<"   Sequences must be in FastA or FastQ format"                    << endl
    << endl
    <<"OUTPUT"                                                           << endl
    <<"   - The output is a tab-separated table containing the columns:" << endl
    <<"     1. id: an internal identifier, provided only for debugging"  << endl
    <<"     2. kmer: the (poly-)nucleotide string"                       << endl
    <<"     3. count: the total count, which can be an approximation"    << endl
    <<"        in scientific notation for large numbers"                 << endl
    <<"   - The extended output is intended to provide higher precision" << endl
    <<"     for large counts, and contains the following additional"     << endl
    <<"     columns (using -e):"                                         << endl
    <<"     4. base: count expressed as times of order of magnitude"     << endl
    <<"     5. order of magnitude: units of base"                        << endl
    <<"     6. extra counts: excess counts from base"                    << endl
    <<"     The actual count can be calculated more accurately as the"   << endl
    <<"     product of base (fourth column) and order of magnitude"      << endl
    <<"     (fourth column), plus the extra counts (fifth column)."      << endl
    <<"     In this compilation, the value of order of magnitude is a"   << endl
    <<"     power of " << OMAG_STEP                                      << endl
    << endl;
   exit(1);
}

void kmer2label(unsigned int &label, char *kmer, unsigned int length){
   // Vars
   char		hex[length+1]; // char representation of some hexadecimal integer

   // label
   for(unsigned int i=0; i<length; i++){
      hex[i] = 'f'; // <- Any other weird thing
      for(int j=0; j<16; j++)
	 if(kmer[i]==nuc_char[j])
	    hex[i]=hex_char[j];
   }
   hex[length] = (char)NULL;
   sscanf(hex, "%x", &label); // hex ---> label
}

void label2kmer(char *&kmer, unsigned int label, unsigned int length){
   // Vars
   char		hex[length+1]; // char representation of some hexadecimal integer

   // label
   sprintf(hex, "%x", label); // hex <--- label
   int		lead_zeroes = length-strlen(hex);
   for(unsigned int i=0; i<lead_zeroes; i++)
      kmer[i] = nuc_char[0];
   for(unsigned int i=lead_zeroes; i<length; i++){
      for(int j=0; j<16; j++)
	 if(hex[i-lead_zeroes]==hex_char[j])
	    kmer[i]=nuc_char[j];
   }
   kmer[length] = (char)NULL;
}

void count_plus(count_t &count, unsigned int howmany){
   // ++hang
   count.hang+=howmany;
   // ++base
   if(count.hang >= count.omag){
      count.base += (unsigned int) count.hang / count.omag;
      count.hang  = count.hang % count.omag;
   }
   // ++omag
   if(count.base >= UINT_MAX/OMAG_STEP){
      count.omag *= OMAG_STEP;
      count.hang += count.base % count.omag;
      count.base /= OMAG_STEP;
   }
}

unsigned int count_polynucleotides(count_t *&counts, char *seqFile, unsigned int k, float heur){
   // Vars
   unsigned int	labels_no, linelen, label;
   ifstream	fileh;
   char		*line, *kmer;

   // Kmers and counts initialization
   labels_no = (unsigned int)pow(16, k);DEBUG(265);
   if(V>=4) fprintf(stderr, "Creating unique integer ids for %u kmers\n", labels_no); 
   counts = (count_t *)malloc(labels_no * (sizeof *counts));
   if(!counts) error("Impossible to allocate memory for that many different kmers", labels_no);
   for(unsigned int a=0; a<labels_no; a++){
      counts[a].kmer = (char *)malloc(k+1);
      if(!counts[a].kmer) error("Insufficient memory to represent the next kmer", a);
      label2kmer(counts[a].kmer, a, k);
      counts[a].label = a;
      counts[a].base = 0;
      counts[a].hang = 0;
      counts[a].omag = 1;
      if(V>=5) fprintf(stderr, "Kmer %s -> %u\n", counts[a].kmer, counts[a].label);
   }
   
   // Allocate char* memory
   if(V>=5) fprintf(stderr, "Allocating %u bits\n", CHAR_BIT*(LARGEST_LINE+1));
   line = (char *)malloc(LARGEST_LINE+1);DEBUG(279);
   if(!line) error("Impossible to allocate one temporal line in memory", LARGEST_LINE);
   if(V>=5) fprintf(stderr, "Allocating %u bits\n", CHAR_BIT*(k+1));
   kmer = (char *)malloc(k+1);
   if(!kmer) error("Impossible to allocate one temporal kmer in memory", k);
   kmer[k] = (char)NULL;

   // Read the file
   if(V>=5) cerr << "Opening file: " << seqFile << endl;
   fileh.open(seqFile, ios::in);DEBUG(286);
   if(!fileh.is_open()) error("Error reading file", seqFile);
   while(!fileh.eof()){
      fileh.getline(line, LARGEST_LINE);
      if(line[0]!='>'){DEBUG(290);
	 linelen = strlen(line);
	 if(linelen < k+1) goto next_line;
	 for(int i=0; i<linelen-k+1; i++){DEBUG(292);
	    memmove(kmer, line+i, k);DEBUG(293);
	    kmer2label(label, kmer, k);DEBUG(294);
	    count_plus(counts[label], 1);DEBUG(295);
	 }
      }
      next_line:;
   }
   fileh.close();

   return labels_no;
}

void report(count_t *&counts, unsigned int length, char *outfile, bool extended){
   ofstream	outfs;
   char		*sep = (char *)"\t", text1[100], text2[100];
   bool		std_out=false;

   if(strlen(outfile)<=0) return;
   else if(strcmp(outfile, "-")==0) std_out = true;
   else{
      outfs.open(outfile, ios::out);
      if(!outfs.is_open()) error("I can not write in the output file", outfile);
   }

   for(unsigned i=0; i<length; i++){
      sprintf(text1, "%u%s%s%s%g",
         counts[i].label, sep,
	 counts[i].kmer, sep,
	 (double)counts[i].base*counts[i].omag+counts[i].hang);
      if(std_out) cout << text1; else outfs << text1;
      if(extended){
         sprintf(text2, "%s%u%s%u%s%u", sep,
	    counts[i].base, sep,
	    counts[i].omag, sep,
	    counts[i].hang);
	 if(std_out) cout << text2; else outfs << text2;
      }
      if(std_out) cout << endl; else outfs << endl;
   }
   if(!std_out) outfs.close();
}

int main(int argc, char *argv[]) {
   cout << "nuc_sampler v1.1" << endl;
   if (argc <= 1) help("");

   // Vars
   char		*file = (char *)"", *format=(char *)"fasta", *outfile=(char *)"-", *namFile, *seqFile;
   double	heur=1.0;
   int		rseed=time(NULL), largest_seq;
   unsigned int	k=1, labels_no, N;
   count_t	*counts;
   bool		extended=false;
   
   // GetOpt
   int		optchr;
   while ((optchr = getopt (argc, argv, "ef:ho:r:s:v:Vk:")) != EOF)
      switch(optchr) {
	 case 'e': extended = true;	break;
	 case 'f': format = optarg;	break;
	 case 'h': help("");		break;
	 case 'k': k = atoi(optarg);	break;
	 case 'o': outfile = optarg;	break;
	 case 'r': rseed=atoi(optarg);	break;
         case 's': file = optarg;	break;
	 case 'v': V = atoi(optarg);	break;
	 case 'V': return 0;
	 case 'x': heur = atof(optarg);	break;
      }
   
   // Initialize
   if (strlen(file) == 0) help("");
   if (strcmp(format, "fasta") != 0 & strcmp(format, "fastq") != 0)
      help("Unsupported value for -f option");
   if (k<=0)
     help("Bad argument for -k option: it must be a non-zero positive integer");
   if (heur<=0.0 || heur>1.0)
     help("Bad argument for -x option: it must be in the range (0, 1]");
   if (outfile && (strlen(outfile)>0) & (strcmp(outfile, "-")!=0))
     remove(outfile);
   srand(rseed);

   // Parse file
   if(V) cerr << "Counting sequences" << endl;
   N = build_index(file, format, namFile, seqFile, largest_seq);
   if(largest_seq<1) error("Your sequences are empty or an internal error occurred.  Largest sequence is ", largest_seq);
   if(V>=2) cerr << "The file " << seqFile << " was just created with sequence-only data" << endl;
   if(V>=4) cerr << "Longest sequence is: " << largest_seq << endl;
   if(N==0) error("The file you provided do not contain sequences.  Before re-run please delete the file", seqFile);
   if(V) cerr << "Reading file with " << N << " sequences" << endl;

   // Run counts
   labels_no = count_polynucleotides(counts, seqFile, k, heur);
   report(counts, labels_no, outfile, extended);

   // Cleanup
   remove(namFile);
   remove(seqFile);

   return 0;
}

