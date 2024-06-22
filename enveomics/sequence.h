
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @license artistic 2.0
// @version 1.0

#ifndef ENVEOMICS_SEQUENCE_H
#define ENVEOMICS_SEQUENCE_H
// #define ENVEOMICS_NUC_T_DEFINE

/**
 * size_t count_seqs(char *file[, const char *format][, int &largest_line[, double &avg_seq]]);
 * Description:
 *   Counts the number of sequences in the file, and optionally measures the maximum and the average length of the
 *   sequences.
 * Input:
 *   char *file: Char array with the path to the file.
 *   char *format (optional): Format of the file.  It can be "fasta", "fastq" or "enveomics-seq".  By default "fasta".
 *   int &largest_line (optional): If passed, saves the length of the largest sequence here.
 *   double &avg_seq (optional): If passed, saves the average sequence length.
 * Output:
 *   size_t: Number of sequences in the dataset.
 */
size_t count_seqs(char *file, const char *format, int &largest_line, double &avg_seq);
size_t count_seqs(char *file, const char *format, int &largest_line);
size_t count_seqs(char *file, const char *format);
size_t count_seqs(char *file, int &largest_line, double &avg_seq);
size_t count_seqs(char *file, int &largest_line);
size_t count_seqs(char *file);

/**
 * size_t build_index(char *sourceFile, char* format, char *&namFileOut, char *&seqFileOut[, int &largest_seq[, double &avg_seq]]);
 * Description:
 *   Builds an index (enveomics-seq format) for the input file.
 * Input:
 *   char *sourceFile: Array of chr with the path to the file.
 *   char *format: Format of the file.  It can be "fasta" or "fastq".
 *   char *&namFileOut: Array of chr to be filled with the path of the output file containing the IDs.
 *   char *&seqFileOut: Array of chr to be filled with the path of the output file containing the sequences.
 *   int &largest_seq (optional): If passed, saves the length of the largest sequence here.
 *   double &avg_seq (optional): If passed, saves the average sequence length here.
 */
size_t build_index(char *sourceFile, char* format, char *&namFileOut, char *&seqFileOut, int &largest_seq, double &avg_seq);
size_t build_index(char *sourceFile, char* format, char *&namFileOut, char *&seqFileOut, int &largest_seq);
size_t build_index(char *sourceFile, char* format, char *&namFileOut, char *&seqFileOut);

/**
 * size_t sub_sample_seqs(char *sourceFile, char *destFile, double portion[, char *format]);
 * Description:
 *   Creates a random sample of sequences (without replacement) at the given portion.
 * Input:
 *   char *sourceFile: Path to the input file.
 *   char *destFile: Path to the output file.
 *   double portion: Portion of the dataset (from 0 to 1) to be sampled.
 *   char *format (optional): Format of the input file.  It can be "fasta", "enveomics-seq" or "fastq".
 *      By default "fasta".
 */
size_t sub_sample_seqs(char *sourceFile, char *destFile, double portion, char *format);
size_t sub_sample_seqs(char *sourceFile, char *destFile, double portion);

/**
 * int get_seqs(char **&seqs, char *file, int from, int number, int largest_seq[, const char *format]);
 * Description:
 *   Takes a file containing sequences and stores in memory starting on the <from>-th sequence as many as
 *   <number> sequences.
 * Input:
 *   char **seq: An array of char arrays to be filled with sequences.
 *   char *file: The file containing the sequences.
 *   int from: First sequence to be saved.
 *   int number: Maximum number of sequences to be saved.
 *   int largest_seq: Length of the largest sequence.
 *   char *format (optional): Format of the file.  It can be "fasta" or "enveomics-seq".  By default "fasta".
 * Output:
 *   Returns the actual number of stored sequences (it can be smaller than <number> if the eof is reached
 *   before.
 */
int get_seqs(char **&seqs, char *file, int from, int number, int largest_seq, char *format);
int get_seqs(char **&seqs, char *file, int from, int number, int largest_seq);

int reverse_complement(char *&out, char *in);
int reverse_complement(string &out, string in);

/*
 * bool has_gz_ext(const char *file);
 * Description:
 *   Evaluates if the input file has a .gz extension
 * Input:
 *   char *file: The path to the input file
 * Output:
 *   Returns true if the file ends in .gz, false otherwise
 */
bool has_gz_ext(const char *file);

/*
 * void unzip_file(char *infile, char *outfile);
 * Description:
 *   Unzips infile it into outfile
 * Input:
 *   char *infile: The path to the input file
 *   char *outfile: The path to the output file
 */
void gunz_file(const char *infile, const char *outfile);

#ifdef ENVEOMICS_NUC_T_DEFINE
#define ENVEOMICS_NUC_T
#include <bitset>
typedef std::bitset<2> nuc_t;
struct nucseq_t {
   nuc_t	*seq;
   size_t	len;
};
nuc_t ctonuc(char c);
char nuctoc(nuc_t nuc);
int atonucseq(nucseq_t &nucseq, char *seq);
int nucseqtoa(char *&charseq, nucseq_t nucseq);
int reverse_complement(nucseq_t &out, nucseq_t in);
#endif

#endif
