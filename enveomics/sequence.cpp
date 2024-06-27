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
#include <sys/types.h>

#include "universal.h"
#include "sequence.h"

using namespace std;

size_t count_seqs(
      char *file, const char *format, int &largest_line, double &avg_seq) {
  // Vars
  unsigned int N = 0, nline = 0, totlen = 0, lastn = 0;
  ifstream fileh;
  int      maxlen = 0;
  char     start;
  bool     isFastQ = false;

  avg_seq = 0.0;

  // Format
  if (strcmp(format, "fasta") == 0) { start = '>'; }
  else if (strcmp(format, "fastq") == 0) { start = '@'; isFastQ=true; }
  else if (strcmp(format, "enveomics-seq") == 0){ start = '>'; }
  else { error("Unsupported format", format); }

  // Count
  fileh.open(file, ios::in);
  if (!fileh.is_open()) error("Impossible to open the file", file);
  while (fileh.good()) {
    string line;
      getline(fileh, line);
      if (line.length() > (size_t) maxlen) maxlen = line.length();
      if ((isFastQ && (nline % 4 == 0)) ||
          (!isFastQ && (line[0] == start))) N++;
      else totlen += line.length();
      if (totlen > UINT_MAX / 10){
        avg_seq = ((avg_seq / N) * lastn) + ((double) totlen / N);
        totlen = 0;
        lastn = N;
      }
      nline++;
  }
  fileh.close();

  largest_line = maxlen;
  avg_seq = ((avg_seq / N) * lastn) + ((double) totlen / N);
  return N;
}

size_t count_seqs(char *file, const char *format, int &largest_line) {
  double dummy;
  return count_seqs(file, format, largest_line, dummy);
}

size_t count_seqs(char *file, const char *format) {
  int    dummy1;
  double dummy2;
  return count_seqs(file, format, dummy1, dummy2);
}

size_t count_seqs(char *file, int &largest_line, double &avg_seq) {
  return count_seqs(file, "fasta", largest_line, avg_seq);
}

size_t count_seqs(char *file, int &largest_line) {
  double dummy;
  return count_seqs(file, largest_line, dummy);
}

size_t count_seqs(char *file) {
  int    dummy1;
  double dummy2;
  return count_seqs(file, dummy1, dummy2);
}

size_t build_index(
      char *sourceFile, char *format, char *&namFileOut, char *&seqFileOut,
      int &largest_seq, double &avg_seq) {
  // Vars
  unsigned int N = 0, nline = 0, totlen = 0, lastn = 0;
  int      maxlen = 0;
  char     start, *namFile, *seqFile;
  bool     inSeq = false, isFastQ = false;
  string   seq, name;
  ifstream infileh, testseq, testnam;
  ofstream seqfileh, namfileh;

  // Format
  if (strcmp(format, "fasta") == 0) {
    start = '>';
  } else if (strcmp(format, "fastq") == 0) {
    start = '@';
    isFastQ = true;
  } else {
    error("Unsupported format", format);
  }

  // Files
  int fileLen = strlen(sourceFile) + 20;
  seqFile = new char[fileLen];
  snprintf(seqFile, fileLen, "%s.enve-seq.%d", sourceFile, getpid());
  namFile = new char[fileLen];
  snprintf(namFile, fileLen, "%s.enve-nam.%d", sourceFile, getpid());
  namFileOut = namFile;
  seqFileOut = seqFile;

  // Check if I need to do it
  testseq.open(seqFile, ios::in);
  if (testseq.is_open()) {
    testseq.close();
    testnam.open(namFile, ios::in);
    if (testnam.is_open()) {
      testnam.close();
      N = count_seqs(seqFile, "enveomics-seq", largest_seq, avg_seq);
      return N;
    }
  }

  // Open file streams
  infileh.open(sourceFile, ios::in);
  if (!infileh.is_open()) error("Cannot open the file", sourceFile);
  seqfileh.open(seqFile, ios::out);
  if (!seqfileh.is_open()) error("Cannot open the file", seqFile);
  namfileh.open(namFile, ios::out);
  if (!namfileh.is_open()) error("Cannot open the file", namFile);

  // Run
  avg_seq = 0.0;
  while (!infileh.eof()) {
    string line;
    getline(infileh, line);
    if ((isFastQ && (nline % 4 == 0)) || (!isFastQ && (line[0] == start))) {
      if (seq.length() > 0) {
        namfileh << ">" << ++N << endl << name << endl;
        seqfileh << ">" <<   N << endl <<  seq << endl;
        if (seq.length() > (size_t) maxlen) maxlen = seq.length();
        totlen += seq.length();
        if (totlen > UINT_MAX/10) {
          avg_seq = ((avg_seq / N) * lastn) + ((double) totlen / N);
          totlen = 0;
          lastn = N;
        }
      }
      name = line.length() > 1 ? line.substr(1) : "";
      seq = (string)"";
      inSeq = true;
    } else {
      if (inSeq) seq.append(line);
      if (isFastQ) inSeq = false;
    }
    nline++;
  }
  if (seq.length() > 0) {
    namfileh << ">" << ++N << endl << name << endl;
    seqfileh << ">" <<   N << endl <<  seq << endl;
    if (seq.length() > (size_t) maxlen) maxlen = seq.length();
    totlen += seq.length();
  }

  // Close file streams
  infileh.close();
  seqfileh.close();
  namfileh.close();

  largest_seq = maxlen;
  avg_seq = ((avg_seq / N) * lastn ) + ((double) totlen / N);
  return N;
}

size_t build_index(
      char *sourceFile, char* format, char *&namFileOut, char *&seqFileOut,
      int &largest_seq) {
  double dummy;
  return build_index(
    sourceFile, format, namFileOut, seqFileOut, largest_seq, dummy
  );
}

size_t build_index(
      char *sourceFile, char* format, char *&namFileOut, char *&seqFileOut) {
  int    dummy1;
  double dummy2;
  return build_index(
    sourceFile, format, namFileOut, seqFileOut, dummy1, dummy2
  );
}

size_t sub_sample_seqs(
      char *sourceFile, char *destFile, double portion, char *format) {
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
      if ((line[0] == start) || !filein.good()) {
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
   for(size_t a=0; a<(size_t)number; a++){
      seqs[a] = new char[largest_seq+1];
      if(!seqs[a]) error("Impossible to allocate memory for another sequence", (unsigned int)a);
   }

   // Read the file
   while(filein.good()){
      string line;
      getline(filein, line);
      if ((line[0] == start) || !filein.good()) {
	 if(entry.size()>0){
	    i++;
	    if(i>=from){
	       if(entry.length() > (size_t)largest_seq) error("Found a sequences largest than expected", (int)entry.length());
	       for(size_t a=0; a<=entry.length(); a++) seqs[n][a] = entry[a];
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

int reverse_complement(string &out, string in){
  int len = in.length();
  out = "";
  for(int i = len-1; i >= 0; i--) {
    out.push_back(in[i]=='A'?'T':
     in[i]=='C'?'G':
     in[i]=='G'?'C':
     in[i]=='T'?'A':
          'N');
  }
  return len;
}

bool has_gz_ext(const char *file) {
  char *ext = new char[4];
  for (int i = 0; i < 4; i++) ext[i] = *(file + strlen(file) - 3 + i);
  return(strcmp(ext, ".gz") == 0);
}

void gunz_file(const char *infile, const char *outfile) {
  char buf[1024 * 1024 + 1];
  int err, len;

  gzFile fi = gzopen(infile, "rb");
  ofstream fo;

  fo.open(outfile);
  gzrewind(fi);
  while (!gzeof(fi)) {
    len = gzread(fi, buf, sizeof(buf) - 1);
    if (len <  0) error(gzerror(fi, &err));
    if (len == 0) error("Unexpected end of file while decompressing");
    buf[len] = '\0';
    fo << buf;
  }
  fo.close();
  gzclose(fi);
}

#ifdef ENVEOMICS_NUC_T
// The bitset representation of nucleotides (in 2 bits)
nuc_t ctonuc(char c) {
  nuc_t nuc;
  nuc.set(); // 11
  if (c == 'A') return (~nuc);     // 00
  if (c == 'C') return (nuc >> 1); // 01
  if (c == 'G') return (nuc << 1); // 10
  if (c == 'T') return nuc;        // 11

  char letter[] = {c, (char)NULL};
  error("Impossible to interpret char as nucleotide", letter);
}

char nuctoc(nuc_t nuc) {
  char c;
  nuc_t nucA = ctonuc('A');
  if (nuc == nucA)        return 'A';
  if (nuc == (nucA >> 1)) return 'C';
  if (nuc == (nucA << 1)) return 'G';
                          return 'T';
}

int atonucseq(nucseq_t &nucseq, char *charseq) {
  int len = strlen(charseq);
  for (size_t a = 0; a < len; a++) nucseq.seq[a] = ctonuc(charseq[a]);
  nucseq.len = len;
  return nucseq.len;
}

int nucseqtoa(char *&charseq, nucseq_t nucseq) {
  for (size_t a = 0; a < nucseq.len; a++) charseq[a] = nuctoc(nucseq.seq[a]);
  charseq[nucseq.len] = (char)NULL;
  return nucseq.len;
}

int reverse_complement(nucseq_t &out, nucseq_t in){ 
  for (size_t i = in.len; i > 0; i--) out.seq[in.len-i] = (~in.seq[i]);
  out.len = in.len;
  return in.len;
}

int seqtoa(char *&charseq, nucseq_t nucseq) {
  for (size_t i = 0; i < nucseq.len; i++) charseq[i] = nuctoc(nucseq.seq[i]);
  charseq[nucseq.len] = (char)NULL;
}

#endif
