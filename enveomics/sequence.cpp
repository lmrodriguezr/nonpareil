// enveomics/sequence.h - Library for sequences on enve-omics software
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @license artistic 2.0
// @version 1.0

#include <iostream>
#include <unistd.h>
#include <fstream>
#include <string>
#include <bitset>
#include <stdio.h>

#include "universal.h"
#include "sequence.h"

using namespace std;

size_t count_seqs(char *file, const char *format, int &largest_line){
   // Vars
   ifstream	fileh;
   unsigned int	N=0;
   int		maxlen=0;
   char		start;
   
   // Format
   if(strcmp(format, "fasta")==0) {start = '>';}
   else if(strcmp(format, "fastq")==0) {start = '@';}
   else if(strcmp(format, "enveomics-seq")==0){start = '>';}
   else {error("Unsupported format", format);}
   
   // Count
   fileh.open(file, ios::in);
   if(!fileh.is_open()) error("Impossible to open the file", file);
   while(fileh.good()){
      string line;
      getline(fileh, line);
      if(line.length() > maxlen) maxlen = line.length();
      if(line[0]==start) N++; 
      if(N==UINT_MAX-1) error("Unable to represent the number of reads, limit reached", UINT_MAX-1);
   }
   fileh.close();
   
   largest_line = maxlen;
   return N;
}

size_t count_seqs(char *file, const char *format){
   int	dummy;
   return count_seqs(file, format, dummy);
}

size_t count_seqs(char *file, int &largest_line){
   return count_seqs(file, "fasta", largest_line);
}

size_t count_seqs(char *file){
   int dummy;
   return count_seqs(file, dummy);
}


size_t build_index(char *sourceFile, char* format, char *&namFileOut, char *&seqFileOut, int &largest_seq){
   // Vars
   unsigned int	N=0;
   int		maxlen=0;
   char		start, end, *namFile, *seqFile;
   bool		inSeq=false;
   string	seq, name;
   ifstream	infileh, testseq, testnam;
   ofstream	seqfileh, namfileh;
   
   // Format
   if(strcmp(format, "fasta")==0) {start = '>'; end='\n';}
   else if(strcmp(format, "fastq")==0) {start = '@'; end='+';}
   else {error("Unsupported format", format);}
   
   // Files
   seqFile = new char[strlen(sourceFile)+10];
   sprintf(seqFile, "%s.enve-seq", sourceFile);
   namFile = new char[strlen(sourceFile)+10];
   sprintf(namFile, "%s.enve-nam", sourceFile);
   namFileOut = namFile;
   seqFileOut = seqFile;

   // Check if I need to do it
   testseq.open(seqFile, ios::in);
   if(testseq.is_open()){
      testseq.close();
      testnam.open(namFile, ios::in);
      if(testnam.is_open()){
	 testnam.close();
         N = count_seqs(seqFile, "enveomics-seq", largest_seq);
	 return N;
      }
   }

   // Open file streams
   infileh.open(sourceFile, ios::in);
   if(!infileh.is_open()) error("Cannot open the file", sourceFile);
   seqfileh.open(seqFile, ios::out);
   if(!seqfileh.is_open()) error("Cannot open the file", seqFile);
   namfileh.open(namFile, ios::out);
   if(!namfileh.is_open()) error("Cannot open the file", namFile);

   // Run
   while(!infileh.eof()){
      string line;
      getline(infileh, line);
      if(line[0]==start){
         if(seq.length()>0){
	    namfileh << ">" << ++N << endl << name << endl;
	    seqfileh << ">" <<   N << endl <<  seq << endl;
	    if(seq.length() > maxlen) maxlen = seq.length();
	    if(N==UINT_MAX-1) error("Impossible to represent the number of reads while indexing, limit reached", UINT_MAX-1);
	 }
         name = line.substr(1);
	 seq = (string)"";
	 inSeq = true;
      }else if(line[0]==end) inSeq = false;
      else if(inSeq) seq.append(line);
   }
   if(seq.length()>0){
      namfileh << ">" << ++N << endl << name << endl;
      seqfileh << ">" <<   N << endl <<  seq << endl;
      if(seq.length() > maxlen) maxlen = seq.length();
   }

   // Close file streams
   infileh.close();
   seqfileh.close();
   namfileh.close();

   largest_seq = maxlen;
   return N;
}

size_t build_index(char *sourceFile, char* format, char *&namFileOut, char *&seqFileOut){
   int dummy;
   return build_index(sourceFile, format, namFileOut, seqFileOut, dummy);
}

size_t sub_sample_seqs(char *sourceFile, char *destFile, double portion, char *format){
   ifstream	filein;
   ofstream	fileout;
   char		start;
   size_t	n=0;
   string	entry;
   
   if(strcmp(format, "fasta")==0) {start = '>';}
   else if(strcmp(format, "enveomics-seq")==0) {start = '>';}
   else if(strcmp(format, "fastq")==0) {start = '@';}
   else {error("Unsupported format", format);}
   
   filein.open(sourceFile, ios::in);
   if(!filein.is_open()) error("Impossible to open the input file", sourceFile);
   fileout.open(destFile, ios::out);
   if(!fileout.is_open()) error("Impossible to open the output file", destFile);
   
   while(filein.good()){
      string line;
      getline(filein, line);
      if(line[0]==start | !filein.good()){
	 if((entry.size()>0) && (portion>=1 || ((double)rand()/RAND_MAX <= portion))){
	    n++;
	    fileout << entry;
	 }
         entry = (string)"";
      }
      entry.append(line);
      entry.append((char *)"\n");
   }

   fileout.close();
   filein.close();
   
   return n;
}

size_t sub_sample_seqs(char *sourceFile, char *destFile, double portion){
   return sub_sample_seqs(sourceFile, destFile, portion, (char *)"fasta");
}

int get_seqs(char **&seqs, char *file, int from, int number, int largest_seq, char *format){
   //  Vars
   ifstream	filein;
   char		start;
   int		n=0, i=0;
   string	entry;

   // Init
   if(strcmp(format, "fasta")==0) {start = '>';}
   else if(strcmp(format, "enveomics-seq")==0) {start = '>';}
   else {error("Unsupported format", format);}
   
   // Open file
   filein.open(file, ios::in);
   if(!filein.is_open()) error("Impossible to open the input file", file);
   
   // Memory allocation
   seqs = new char*[number];
   if(!seqs) error("Impossible to allocate memory for that many sequences", number);
   for(unsigned int a=0; a<number; a++){
      seqs[a] = new char[largest_seq+1];
      if(!seqs[a]) error("Impossible to allocate memory for another sequence", a);
   }
   
   // Read the file
   while(filein.good()){
      string line;
      getline(filein, line);
      if(line[0]==start | !filein.good()){
	 if(entry.size()>0){
	    i++;
	    if(i>=from){
	       if(entry.length() > largest_seq) error("Found a sequences largest than expected", (int)entry.length());
	       for(int a=0; a<=entry.length(); a++) seqs[n][a] = entry[a];
	       n++;
	       if(n >= number) break;
	    }
	 }
	 entry = (string)"";
      }else{
         entry.append(line);
         //entry.append((char *)"\n");
      }
   }
   
   // Finalize
   filein.close();
   return n;
}

int get_seqs(char **&seqs, char *file, int from, int number, int largest_seq){
   return get_seqs(seqs, file, from, number, largest_seq, (char *)"fasta");
}

int reverse_complement(char *&out, char *in){
   int	len = strlen(in);
   for(int i=len; i>0; i--)
      out[len-i] = in[i]=='A'?'T':
		   in[i]=='C'?'G':
		   in[i]=='G'?'C':
		   in[i]=='T'?'A':
			      'N';
   out[len]=(char)NULL;
   return len;
}

#ifdef ENVEOMICS_NUC_T
// The bitset representation of nucleotides (in 2 bits)
nuc_t ctonuc(char c){
   nuc_t nuc;
   nuc.set(); // 11
   if(c=='A') return (~nuc);	// 00
   if(c=='C') return (nuc>>1);	// 01
   if(c=='G') return (nuc<<1);	// 10
   if(c=='T') return nuc;	// 11

   char	letter[] = {c, (char)NULL};
   error("Impossible to interpret char as nucleotide", letter);
}

char nuctoc(nuc_t nuc){
   char c;
   nuc_t nucA = ctonuc('A');
   if(nuc==nucA)	return 'A';
   if(nuc==(nucA>>1))	return 'C';
   if(nuc==(nucA<<1))	return 'G';
   			return 'T';
}

int atonucseq(nucseq_t &nucseq, char *charseq){
   int	len = strlen(charseq);
   for(size_t a=0; a<len; a++) nucseq.seq[a] = ctonuc(charseq[a]);
   nucseq.len=len;
   return nucseq.len;
}

int nucseqtoa(char *&charseq, nucseq_t nucseq){
   for(size_t a=0; a<nucseq.len; a++) charseq[a] = nuctoc(nucseq.seq[a]);
   charseq[nucseq.len] = (char)NULL;
   return nucseq.len;
}

int reverse_complement(nucseq_t &out, nucseq_t in){
   for(size_t i=in.len; i>0; i--) out.seq[in.len-i] = (~in.seq[i]);
   out.len = in.len;
   return in.len;
}

int seqtoa(char *&charseq, nucseq_t nucseq){
   for(size_t i=0; i<nucseq.len; i++) charseq[i] = nuctoc(nucseq.seq[i]);
   charseq[nucseq.len] = (char)NULL;
}

#endif

