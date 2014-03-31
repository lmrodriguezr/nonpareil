// nonpareil_mating - Part of the nonpareil package
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @license artistic 2.0
// @version 1.0

#ifndef ENVEOMICS_NONPAREIL_MATING_H
#define ENVEOMICS_NONPAREIL_MATING_H

#include <fstream>

/**
 * typedef matepar_t;
 * Description:
 *   Parameters of matting.
 * Elements:
 *   double overlap: Minimum overlap of the two reads, as a fraction of the shortest.
 *   double similarity: Minimum similarity, as a fraction of the overlapping region.
 *   double qryportion: Portion of the dataset to be used as query set.
 *   bool revcom: True if the sequences must be matted in both strands.
 *   bool n_as_mismatch: True if the 'N' characters must be treated as a mismatch.
 */
typedef struct {
   double	overlap;
   double	similarity;
   double	qryportion;
   bool		revcom;
   bool		n_as_mismatch;
} matepar_t;

/**
 * typedef matejob_t;
 * Description:
 *   Arguments required by on threaded matting job.
 * Elements:
 *   int id: The consecutive number of the thread.
 *   size_t from: Initial element of the block to mate (start of the region).
 *   size_t number: Number of elements in the block to mate (size of the region).
 *   int from_in_result: Position of the initial element of the block in the query set (start of the block).
 *   matepar_t par: Parameters of the matting.
 *   int **results: Reference to an int array to be filled with the matting results.
 *   pthread_mutex_t *mutex: Reference to a mutex to be used when modifying the array pointed by result.
 *   char ***blockA: Reference to an array of char arrays containing the reads in Block A.
 *   char ***blockB: Reference to an array of char arrays containing the reads in Block B.
 *   int size_blockA: Number of elements in the array pointed by blockA.
 *   int size_blockB: Number of elements in the array pointed by blockB.
 */
typedef struct {
   // For the thread
   int			id;
   size_t 		from;
   size_t 		number;
   // Universal
   int			from_in_result;
   matepar_t		par;
   int			**result;
   pthread_mutex_t	*mutex;
   // The job
   char			***blockA;
   char			***blockB;
   int			size_blockA;
   int			size_blockB;
} matejob_t;

/**
 * size_t nonpareil_mate(int *&result, char *file[, char *q_file]
 *			int threads, unsigned int lines_in_ram,
 *			unsigned int total_seqs, unsigned int largest_seq,
 *			[unsigned int q_largest_seq,]
 *			matepar_t matepar);
 * Description:
 *   Counts all the mates in the dataset by blocks determined by the ammount of RAM allowed to use.
 * Input:
 *   int *&result: Array of integers to be filled with the results of mating.
 *   char *file: File containing the sequences (assumed to be in "enveomics-seq" format).
 *   char *q_file (optional): File containing the query sequences (assumed to be in "enveomics-seq" format).
 *		If not passed, `file` is used as query file.
 *   int threads: Number of threads (CPUs) to use.
 *   unsigned int lines_in_ram: Maximum number of sequences to be stored in RAM.
 *   unsigned int total_seqs: Total number of sequences in the file.
 *   unsigned int largest_seq: Length of the largest sequence.
 *   unsigned int q_largest_seq (optional): Length of the largest query sequence. Required if `q_file` is
 		passed.
 *   matepar_t matepar: Parameters for the mating.
 * Output:
 *   size_t: The number of slots in the results array (equal to the number of sampled queries).
 */
size_t nonpareil_mate(int *&result, char *file,
			int threads, unsigned int lines_in_ram,
			unsigned int total_seqs, unsigned int largest_seq,
			matepar_t matepar);
size_t nonpareil_mate(int *&result, char *file, char *q_file,
			int threads, unsigned int lines_in_ram,
			unsigned int total_seqs, unsigned int largest_seq,
			unsigned int q_largest_seq,
			matepar_t matepar);

/**
 * void nonpareil_count_mates_block(int *&result, int from_in_result,
 *			char **&blockA, char **&blockB, int sizeBlockA, int sizeBlockB,
 *			int threads, matepar_t matepar);
 * Description:
 *   Performs the matting of the reads in one single block.
 * Input:
 *   int *&result: An int array to be filled with the number of mates per query read.
 *   int from_in_result: Position of the initial read in blockA, with respect to
 *      all the reads in the query set.  Also, the initial position in the result array
 *      to be filled.
 *   char **&blockA: Array of chr arrays containing the block of query sequences to mate.
 *   char **&blockB: Array of chr arrays containing the block of subject sequences to mate.
 *   int sizeBlockA: Number of elements in blockA.
 *   int sizeBlockB: Number of elements in blockB.
 *   int threads: Number of threads to use.
 *   matepar_t matepar: Parameters of the matting passed in a matepar_t complete object.
 */
void nonpareil_count_mates_block(int *&result, int from_in_result,
			char **&blockA, char **&blockB, int sizeBlockA, int sizeBlockB,
			int threads, matepar_t matepar);

/**
 * void *nonpareil_count_mates_thr(void *matejob_ref);
 * Description:
 *   Interface to nonpareil_count_mates() intended to allow threaded search.
 * Input:
 *   void *matejob_ref: Reference to a complete matejob_t object.
 * Output:
 *   void *: Reference to a zero.  Otherwise, reference to an error message.
 */
void *nonpareil_count_mates_thr(void *matejob_ref);

/**
 * void nonpareil_count_mates(int *&result, char **&blockA, char **&blockB,
 *			int fromA, int numberA, int fromB, int numberB,
 *			int talk, matepar_t matepar);
 * Description:
 *   Performs the matting on a region of a block of sequences.
 * Input:
 *   int *&result: Int array to be filled with the matting results.
 *   char **&blockA: Array of chr arrays containing the block of query sequences to mate.
 *   char **&blockB: Array of chr arrays containing the block of subject sequences to mate.
 *   int fromA: First element of blockA to be matted.
 *   int numberA: Number of elements in blockA to be matted.
 *   int fromB: First element of blockB to be matted.
 *   int numberB: Number of elements in blockB to be matted.
 *   int talk: The period (in number of query sequences) to report advance via say().
 *      Zero (0) to be silent.
 *   matepar_t matepar: A complete matepar_t object containing the matting parameters.
 */
void nonpareil_count_mates(int *&result, char **&blockA, char **&blockB,
			int fromA, int numberA, int fromB, int numberB,
			int talk, matepar_t matepar);

/**
 * bool nonpareil_compare_reads(char *seqA, char *seqB, matepar_t matepar);
 * Description:
 *   Compares two reads.
 * Input:
 *   char *seqA: Read A.
 *   char *seqB: Read B.
 *   matepar_t matepar: Parameters of the comparison.
 * Output:
 *   bool: True if the two sequences are similar, false otherwise.
 */
bool nonpareil_compare_reads(char *seqA, char *seqB, matepar_t matepar);

/**
 * bool nonpareil_compare_reads_shortfirst(char *seqA, char *seqB, int lenA, int lenB, matepar_t matepar);
 * Description:
 *   Compares two reads, the first of which is shorter or equal in length to the second.
 * Input:
 *   char *seqA: Read A.
 *   char *seqB: Read B.
 *   int lenA: Length of seqA.
 *   int lenB: Length of seqB.  Notice that you MUST ensure that lenA<=lenB, the function doesn't check!
 *   matepar_t matepar: Parameters of the comparison.
 * Output:
 *   bool: True if the two sequences are similar, false otherwise.
 */
bool nonpareil_compare_reads_shortfirst(char *seqA, char *seqB, int lenA, int lenB, matepar_t matepar);

/**
 * void nonpareil_save_mates(int *&result, int no_results, char *file);
 * Description:
 *   Saves a Mates vector into a file.
 * Input:
 *   int *&result: Mates vector.
 *   int no_results: Number of entries in Mates.
 *   char *file: Path to the file where the vector must be saved.
 */
void nonpareil_save_mates(int *&result, int no_results, char *file);

#endif

