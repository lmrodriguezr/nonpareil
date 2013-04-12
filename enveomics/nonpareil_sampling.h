// nonpareil_sampling - Part of the nonpareil package
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @license artistic 2.0
// @version 1.0

#ifndef ENVEOMICS_NONPAREIL_SAMPLING_H
#define ENVEOMICS_NONPAREIL_SAMPLING_H

/**
 * typedef sample_t;
 * Description:
 *   An object containing the summary of a sub-sample.
 * Elements:
 *   double portion: Portion at which the sub-sample was generated.
 *   double avg: Average of the non-unique portion.
 *   double sd: Standard deviation of the non-unique portion.
 *   double q1: Quartile I of the non-unique portion.
 *   double q2: Quartile II of the non-unique portion.
 *   double q3: Quartile III of the non-unique portion.
 */
typedef struct {
   double	portion;
   double	avg;
   double	sd;
   double	q1;
   double	q2;
   double	q3;
} sample_t;

/**
 * typedef samplepar_t;
 * Description:
 *   Collection of parameters for sub-sampling.
 * Elements:
 *   double portion: Portion at which the sub-sample must be generated.
 *   double portion_min: Minimum portion in the run.
 *   double portion_max: Maximum portion in the run.
 *   double portion_itv: Interval of portions in the run.
 *   int replicates: Number of replicates to generate.
 *   int **mates: Reference to the int array containing the matting results.
 *   int mates_size: Number of elements in the array referenced by mates.
 *   int total_reads: Total number of reads in the dataset.
 *   bool portion_as_label: True if the label column of the summary file must
 *      be filled with the portion, instead of the number of reads.
 */
typedef struct {
   double	portion;
   double	portion_min;
   double	portion_max;
   double	portion_itv;
   int		replicates;
   int		**mates;
   int		mates_size;
   int		total_reads;
   bool		portion_as_label;
} samplepar_t;

/**
 * typedef samplejob_t;
 * Description:
 *   An object defining the job to run by a single sample thread.
 * Elements:
 *   int id: The consecutive number of the thread.
 *   int from_in_result: The initial index on the result array to be filled.
 *   int number: The number of sub-samples to perform.
 *   samplepar_t samplepar: The parameters of the sub-sampling.
 *   double **result: A reference to a double array to be filled with the
 *      sub-sampling results.
 *   pthread_mutex_t *mutex: A mutex to be used when modifying the result array.
 */
typedef struct {
   int			id;
   int			from_in_result;
   int			number;
   samplepar_t		samplepar;
   double		**result;
   pthread_mutex_t	*mutex;
} samplejob_t;

/**
 * int nonpareil_sample_portion(double *&result, int threads, samplepar_t samplepar);
 * Description:
 *   Creates an array of sub-samples given a matted set.
 * Input:
 *   double *&result: A pointer to a double array, to be filled with the results of the sub-samples.
 *   int threads: Number of threads to use in the sub-sampling process.
 *   samplepar_t samplepar: Parameters of the sub-sampling.  It's expected to be fully defined.
 * Output:
 *   int: Number of samples taken.  This should be equal to the replicates element of samplepar.
 */
int nonpareil_sample_portion(double *&result, int threads, samplepar_t samplepar);

/**
 * void *nonpareil_sample_portion_thr(void *samplejob_ref);
 * Description:
 *   Performs the sub-sampling in a region of the matting resuts array.  All the parameters are passed
 *   though a reference to a samplejob_t element to allow threaded execution.
 * Input:
 *   void *samplejob_ref: A reference to a complete samplejob_t object.
 * Output:
 *   void *: A reference to a zero.  Otherwise, a reference to an error message.
 */
void *nonpareil_sample_portion_thr(void *samplejob_ref);

/**
 * sample_t nonpareil_sample_summary(double *&sample_result, int sample_number, char *alldata, char *outfile, samplepar_t samplepar);
 * Description:
 *   Calculates and saves the summary of a sub-sampling.
 * Input:
 *   double *&sample_result: A double array containing sub-sampling results, as produced by nonpareil_sample_portion().
 *   int sample_number: Number of elements in sample_result.
 *   char *alldata: A char array with the path to the file in which all-data report must be created.  It can be an empty string (""),
 *      to ignore the all-data report.
 *   char *outfile: Char array with the path to the file in which the summary report must be created.  It can be a dash ("-") to send
 *      the summary to the stdout, or an empty string ("") to ignore the report.
 *   samplepar_t samplepar: A samplepar_t object with the sub-sampling parameters.  It's expected to be complete.
 * Output:
 *   sample_t: A sample_t object containing the summary of the sub-sample.
 */
sample_t nonpareil_sample_summary(double *&sample_result, int sample_number, char *alldata, char *outfile, samplepar_t samplepar);

#endif

