
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @license artistic 2.0
// @version 1.0

#ifndef ENVEOMICS_SEQUENCE_H
#define ENVEOMICS_SEQUENCE_H

#include <vector>

/**
 * size_t count_seqs(
 *       char *file[, const char *format]
 *       [, int &largest_line[, double &avg_seq]]);
 * Description:
 *   Counts the number of sequences in the file, and optionally measures the
 *   maximum and the average length of the sequences.
 * Input:
 *   - `char *file`: Char array with the path to the file
 *   - `char *format` (optional): Format of the file.  It can be "fasta",
 *     "fastq" or "enveomics-seq".  By default "fasta"
 *   - `int &largest_line` (optional): If passed, saves the length of the
 *     largest sequence here
 *   - `double &avg_seq` (optional): If passed, saves the average sequence
 *     length
 * Output:
 *   `size_t`: Number of sequences in the dataset.
 */
size_t count_seqs(
      char *file, const char *format, int &largest_line, double &avg_seq);
size_t count_seqs(char *file, const char *format, int &largest_line);
size_t count_seqs(char *file, const char *format);
size_t count_seqs(char *file, int &largest_line, double &avg_seq);
size_t count_seqs(char *file, int &largest_line);
size_t count_seqs(char *file);

/**
 * size_t build_index(
 *       char *sourceFile, char* format, char *&namFileOut, char *&seqFileOut
 *       [, int &largest_seq[, double &avg_seq[, int len_min[, bool do_nam]]]]);
 * Description:
 *   Builds an index (enveomics-seq format) for the input file.
 * Input:
 *   - `char *sourceFile`: Array of chr with the path to the file.
 *   - `char *format: Format of the file.  It can be "fasta" or "fastq".
 *   - `char *&namFileOut`: Array of chr to be filled with the path of the
 *     output file containing the IDs.
 *   - `char *&seqFileOut`: Array of chr to be filled with the path of the
 *     output file containing the sequences.
 *   - `int &largest_seq` (optional): If passed, saves the length of the largest
 *     sequence here.
 *   - `double &avg_seq` (optional): If passed, saves the average sequence
 *     length here
 *   - `int len_min` (optional): If passed, minimum read length, by default: 1.
 *   - `bool do_nam` (optional): If true (default), saves the original sequence
 *     names in a dedicated fasta-like file. If `false`, the file will still be
 *     created and the path set on `namFileOut`, but it will be an empty file
 */
size_t build_index(
      char *sourceFile, char* format, char *&namFileOut, char *&seqFileOut,
      int &largest_seq, double &avg_seq, int len_min, bool do_nam);
size_t build_index(
      char *sourceFile, char* format, char *&namFileOut, char *&seqFileOut,
      int &largest_seq, double &avg_seq, int len_min);
size_t build_index(
      char *sourceFile, char* format, char *&namFileOut, char *&seqFileOut,
      int &largest_seq, double &avg_seq);
size_t build_index(
      char *sourceFile, char* format, char *&namFileOut, char *&seqFileOut,
      int &largest_seq);
size_t build_index(
      char *sourceFile, char* format, char *&namFileOut, char *&seqFileOut);

/**
 * size_t sub_sample_seqs(
 *       char *sourceFile, char *destFile, double portion[, char *format]);
 * Description:
 *   Creates a random sample of sequences (without replacement) at the given
 *   portion.
 * Input:
 *   - `char *sourceFile`: Path to the input file.
 *   - `char *destFile`: Path to the output file.
 *   - `double portion`: Portion of the dataset (from 0 to 1) to be sampled.
 *   - `char *format` (optional): Format of the input file.  It can be "fasta",
 *     "enveomics-seq" or "fastq". By default "fasta".
 */
size_t sub_sample_seqs(
      char *sourceFile, char *destFile, double portion, char *format);
size_t sub_sample_seqs(
      char *sourceFile, char *destFile, double portion);

/**
 * int get_seqs(
 *       char **&seqs, char *file, int from, int number,
 *       int largest_seq[, const char *format]);
 * Description:
 *   Takes a file containing sequences and stores in memory starting on the
 *   `from`-th sequence as many as <number> sequences.
 * Input:
 *   - `char **seq`: An array of char arrays to be filled with sequences.
 *   - `char *file`: The file containing the sequences.
 *   - `int from`: First sequence to be saved.
 *   - `int number`: Maximum number of sequences to be saved.
 *   - `int largest_seq`: Length of the largest sequence.
 *   - `char *format` (optional): Format of the file.  It can be "fasta" or
 *     "enveomics-seq".  By default "fasta".
 * Output:
 *   Returns the actual number of stored sequences; it can be smaller than
 *   `number` if the EOF is reached before.
 */
int get_seqs(
      char **&seqs, char *file, int from, int number, int largest_seq,
      char *format);
int get_seqs(
      char **&seqs, char *file, int from, int number, int largest_seq);

/**
 * int reverse_complement(char *&out,  char  *in);
 * int reverse_complement(string &out, string in);
 * Description:
 *   Calculate the reverse-complement of a nucleotide sequence
 * Input:
 *   - [`char *`|`string`] `&out`: The output, reverse-complemented sequence
 *   - [`char *`|`string`] `in`: The input sequence
 * Output:
 *   Returns the length of the sequence
 */
int reverse_complement(char *&out, char *in);
int reverse_complement(string &out, string in);

/**
 * std::vector<unsigned int> get_seq_lengths(char *file, unsigned int total_seqs);
 * Description:
 *   Performs a single lightweight pass over a file in "enveomics-seq" format,
 *   returning the length of every sequence, in file order, without loading
 *   the full sequence content into memory.
 * Input:
 *   - `char *file`: The file containing the sequences (in "enveomics-seq"
 *     format).
 *   - `unsigned int total_seqs`: Total number of sequences in the file, used
 *     to reserve the output vector.
 * Output:
 *   Returns a vector with the length of each sequence, in file order.
 */
std::vector<unsigned int> get_seq_lengths(char *file, unsigned int total_seqs);

/**
 * void write_seq_range_to_fasta(
 *       char *file, char *outfile, size_t start, unsigned int count);
 * Description:
 *   Writes sequences `[start, start + count)` (1-based) from a file in
 *   "enveomics-seq" format into a fresh FASTA file.
 * Input:
 *   - `char *file`: The source file (in "enveomics-seq" format).
 *   - `char *outfile`: Path to the FASTA file to create.
 *   - `size_t start`: First sequence to write (1-based).
 *   - `unsigned int count`: Maximum number of sequences to write.
 */
void write_seq_range_to_fasta(
      char *file, char *outfile, size_t start, unsigned int count);

/**
 * bool has_gz_ext(const char *file);
 * Description:
 *   Evaluates if the input file has a .gz extension
 * Input:
 *   - `char *file`: The path to the input file
 * Output:
 *   Returns true if the file ends in .gz, false otherwise
 */
bool has_gz_ext(const char *file);

/**
 * void gunz_file(const char *infile, const char *outfile);
 * Description:
 *   Decompresses a gzip-compressed file in full, writing the result to a
 *   new plain file. Most of this codebase reads `.gz` input directly
 *   (see `build_index`), without needing this; it's still used where a
 *   consumer genuinely needs a fully decompressed copy on disk.
 * Input:
 *   - `char *infile`: The path to the input (gzip-compressed) file
 *   - `char *outfile`: The path to the output (decompressed) file
 */
void gunz_file(const char *infile, const char *outfile);

// The following block (a 2-bit-per-base nucleotide encoding, distinct from
// the "enveomics-seq" file format used elsewhere in this file) is opt-in
// via ENVEOMICS_NUC_T_DEFINE, which nothing in this codebase currently
// defines -- no functions here are compiled or used today. Kept available
// for callers that want a compact in-memory nucleotide representation.
#ifdef ENVEOMICS_NUC_T_DEFINE
#define ENVEOMICS_NUC_T
#include <bitset>

// A single nucleotide, 2 bits: 00=A, 01=C, 10=G, 11=T (see `ctonuc`).
typedef std::bitset<2> nuc_t;

// A nucleotide sequence stored in the 2-bit `nuc_t` encoding.
struct nucseq_t {
   nuc_t	*seq;
   size_t	len;
};

/**
 * nuc_t ctonuc(char c);
 * Description:
 *   Encodes a single nucleotide character ('A'/'C'/'G'/'T') as a `nuc_t`.
 * Input:
 *   - `char c`: The nucleotide character to encode
 * Output:
 *   Returns the corresponding `nuc_t`. Behavior is undefined for
 *   characters other than 'A', 'C', 'G', or 'T'
 */
nuc_t ctonuc(char c);

/**
 * char nuctoc(nuc_t nuc);
 * Description:
 *   Decodes a `nuc_t` back into its nucleotide character.
 * Input:
 *   - `nuc_t nuc`: The encoded nucleotide
 * Output:
 *   Returns 'A', 'C', 'G', or 'T'
 */
char nuctoc(nuc_t nuc);

/**
 * int atonucseq(nucseq_t &nucseq, char *seq);
 * Description:
 *   Encodes a nucleotide-character sequence into a `nucseq_t`. `nucseq.seq`
 *   must already point to a buffer with room for `strlen(seq)` `nuc_t`
 *   elements.
 * Input:
 *   - `nucseq_t &nucseq`: The output, encoded sequence
 *   - `char *seq`: The input sequence, as nucleotide characters
 * Output:
 *   Returns the length of the sequence (same as `nucseq.len` on return)
 */
int atonucseq(nucseq_t &nucseq, char *seq);

/**
 * int nucseqtoa(char *&charseq, nucseq_t nucseq);
 * Description:
 *   Decodes a `nucseq_t` back into a nucleotide-character string.
 *   `charseq` must already point to a buffer with room for `nucseq.len + 1`
 *   characters (including the terminating null).
 * Input:
 *   - `char *&charseq`: The output, decoded sequence (null-terminated)
 *   - `nucseq_t nucseq`: The input, encoded sequence
 * Output:
 *   Returns the length of the sequence (excluding the terminating null)
 */
int nucseqtoa(char *&charseq, nucseq_t nucseq);

/**
 * int reverse_complement(nucseq_t &out, nucseq_t in);
 * Description:
 *   Calculate the reverse-complement of a `nucseq_t`-encoded sequence.
 *   `out.seq` must already point to a buffer with room for `in.len`
 *   `nuc_t` elements.
 * Input:
 *   - `nucseq_t &out`: The output, reverse-complemented sequence
 *   - `nucseq_t in`: The input sequence
 * Output:
 *   Returns the length of the sequence
 */
int reverse_complement(nucseq_t &out, nucseq_t in);
#endif

#endif
