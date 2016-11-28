// nonpareil_mating - Part of the nonpareil package
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @license artistic 2.0
// @version 1.0

#define _MULTI_THREADED
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <sys/types.h>

#include "universal.h"
#include "multinode.h"
#include "sequence.h"
#include "nonpareil_mating.h"

extern int processID;
extern int processes;

#define LARGEST_PATH 4096

using namespace std;

size_t nonpareil_mate(int *&result, char *file,
		int threads, unsigned int lines_in_ram,
		unsigned int total_seqs, unsigned int largest_seq,
		matepar_t matepar){
	// Use `file` for both query and subject sequences
	return nonpareil_mate(result, file, file, threads, lines_in_ram,
			total_seqs, largest_seq, largest_seq, matepar);
}
size_t nonpareil_mate(int *&result, char *file, char *q_file,
		int threads, unsigned int lines_in_ram,
		unsigned int total_seqs, unsigned int largest_seq,
		unsigned int q_largest_seq,
		matepar_t matepar){
   // Vars
   int		no_blocks_sbj=0, no_blocks_qry=0,
   		no_seqs_block_qry=0, no_seqs_block_sbj=0,
		tmp_ram, result_i=0,
		size_blockA, size_blockB;
   size_t	qry_seqs=0;
   char		**blockA, **blockB, *sampleFile;

   // Set subsampling
   //sampleFile = (char *)malloc(LARGEST_PATH * (sizeof *sampleFile));
   sampleFile = new char[LARGEST_PATH];
   if(processID==0){
      sprintf(sampleFile, "%s.subsample.%d", q_file, getpid());
      say("3ss$", "Building query set at ", sampleFile);
      qry_seqs = sub_sample_seqs(q_file, sampleFile, matepar.qryportion, (char *)"enveomics-seq");
      say("4sus$", "Query set built with ", qry_seqs, " sequences");
      if(qry_seqs==0) error("Impossible to create the query set.  Is the -X/-x value too small?");
   }
   sampleFile = broadcast_char(sampleFile, LARGEST_PATH);
   qry_seqs = broadcast_int(qry_seqs);

   // Blank results
   result = new int[qry_seqs];
   for(size_t a=0; a<qry_seqs; a++) result[a] = 0;

   // Design blocks
   if(processID==0){
      say("5sis$", "Designing the blocks scheme for ", total_seqs, " sequences");

      no_blocks_qry = (int)ceil((double)qry_seqs*2/(double)lines_in_ram); // <- Maximum half of the available slots
      if(no_blocks_qry==0) no_blocks_qry=1; // <-- Because of float precision
      no_seqs_block_qry = (int)ceil((double)qry_seqs/(double)no_blocks_qry);
      say("6sisi$", "Qry blocks:", no_blocks_qry, ", seqs/block:", no_seqs_block_qry);

      no_blocks_sbj = (int)ceil( (double)total_seqs/(double)(lines_in_ram - no_seqs_block_qry) );
      if(no_blocks_sbj==0) no_blocks_sbj=1; // <-- Because of float precision
      no_blocks_sbj = (int)ceil( (double)no_blocks_sbj/(double)processes )*processes;
      no_seqs_block_sbj = (int)ceil((double)total_seqs/(double)no_blocks_sbj);
      say("6sisi$", "Sbj blocks:", no_blocks_sbj, ", seqs/block:", no_seqs_block_sbj);
   }
   no_blocks_qry = broadcast_int(no_blocks_qry);
   no_seqs_block_qry = broadcast_int(no_seqs_block_qry);
   no_blocks_sbj = broadcast_int(no_blocks_sbj);
   no_seqs_block_sbj = broadcast_int(no_seqs_block_sbj);

   // Mating
   if(processID==0) say("3sisis$", "Mating sequences in ", no_blocks_qry, " by ", no_blocks_sbj, " blocks");
   if(processID==0 && processes>1) say("3sis$", "Silencing log in slave processes (", processes-1, ")");
   for(int i=0; i<no_blocks_qry; i++){
      // Sequences in block A (qry)
      tmp_ram = (int)(((double)no_seqs_block_qry/1024)*q_largest_seq*(sizeof **blockA)/1024);
      if(processID==0) say("5sisi$", "Allocating ~", tmp_ram, " Mib in RAM for block qry:", i+1);
      size_blockA = get_seqs(blockA, sampleFile, i*no_seqs_block_qry+1, no_seqs_block_qry, q_largest_seq, (char *)"enveomics-seq");
      if(size_blockA==0) error("Impossible to get the i-th query block", i);

      // Sequences in block B (sbj)
      for(int j=0; j<no_blocks_sbj; j++){
	 if(j%processes == processID){
	    tmp_ram = (int)(((double)no_seqs_block_sbj/1024)*largest_seq*(sizeof **blockB)/1024);
	    if(processID==0) say("5sisi$", "Allocating ~", tmp_ram, " Mib in RAM for block sbj:", j+1);
	    size_blockB = get_seqs(blockB, file, j*no_seqs_block_sbj+1, no_seqs_block_sbj, largest_seq, (char *)"enveomics-seq");
	    if(size_blockB==0) error("Impossible to get the i-th subject block", j);

	    // Mate
	    if(processID==0) say("4sisi$", "Computing block ", (i+1)*(j+1), "/", no_blocks_qry*no_blocks_sbj);
	    nonpareil_count_mates_block(result, result_i, blockA, blockB, size_blockA, size_blockB, threads, matepar);
	    for(int a=0; a<size_blockB; a++) delete [] blockB[a];
	    delete[] blockB;
	 }
      }
      result_i += size_blockA;
      for(int a=0; a<size_blockA; a++) delete [] blockA[a];
      delete[] blockA;
   }
   barrier_multinode();

   // Reduce multi-node results
   if(processes>1){
      // DEBUG
      //char *cntfile = new char[123];
      //sprintf(cntfile, "T.c%d", processID);
      //nonpareil_save_mates(result, qry_seqs, cntfile);
      // END DEBUG
      int *result_sum = new int[qry_seqs];
      reduce_sum_int(result, result_sum, qry_seqs);
      if(processID==0) for(size_t a=0; a<qry_seqs; a++) result[a] = result_sum[a];
   }
   barrier_multinode();

   if(processID==0) remove(sampleFile);

   return qry_seqs;
}

void nonpareil_count_mates_block(int *&result, int from_in_result,
			char **&blockA, char **&blockB, int sizeBlockA, int sizeBlockB,
			int threads, matepar_t matepar){
   // Vars
   if(sizeBlockA<threads) threads = sizeBlockA;
   pthread_t		thread[threads];
   matejob_t		matejob[threads];
   pthread_mutex_t	mutex = PTHREAD_MUTEX_INITIALIZER;
   unsigned int		mates_per_thr = (unsigned int)ceil((double)sizeBlockA/threads);
   int			rc;

   // Set threads
   threads = (int)ceil((double)sizeBlockA/mates_per_thr);

   // Launch jobs
   if(processID==0) say("5sis$", "Launching parallel comparisons to ", threads, " threads");
   for(int thr=0; thr<threads; thr++){
     matejob[thr].id = thr; // The ID of the thread
     matejob[thr].from = mates_per_thr*thr; // The first qry sequence to process (in zero-count)
     matejob[thr].number = (matejob[thr].from+mates_per_thr > (size_t)sizeBlockA ? sizeBlockA-matejob[thr].from : mates_per_thr); // How many qry sequences to process
     matejob[thr].from_in_result = from_in_result; // Where to start saving results (in zero-count)
     matejob[thr].par = matepar; // It's cheap to create multiple copies of this, and it's safer than passing a reference.
     matejob[thr].result = &result; // This is only used with mutex (to save RAM)
     matejob[thr].mutex = &mutex; // And this is the mutex, it MUST be the same
     matejob[thr].blockA = &blockA;
     matejob[thr].blockB = &blockB;
     matejob[thr].size_blockA = sizeBlockA;
     matejob[thr].size_blockB = sizeBlockB;

     if(matejob[thr].number==0)
        error("Unexpectedly, the thread contains zero load", thr);

     if((rc=pthread_create(&thread[thr], NULL, &nonpareil_count_mates_thr, (void *)&matejob[thr]  )))
        error("Thread creation failed", (char)rc);
   }

   // Gather jobs
   for(int thr=0; thr<threads; thr++){
      pthread_join(thread[thr], NULL);
   }

   return;
}

void *nonpareil_count_mates_thr(void *matejob_ref){
   // Vars
   matejob_t	*matejob = (matejob_t *)matejob_ref;
   int		*result_cp = new int[matejob->number];
   if(!result_cp) error("Impossible to allocate memory for the results of the new thread", matejob->id);

   // Run comparisons
   for(size_t a=0; a<matejob->number; a++) result_cp[a] = 0;
   nonpareil_count_mates(result_cp, *matejob->blockA, *matejob->blockB,
   	matejob->from, matejob->number, 0, matejob->size_blockB,
	(matejob->id==0?(int)ceil((double)matejob->number/100.0):0), matejob->par);

   // Transfer the results to the external array
   pthread_mutex_lock( matejob->mutex );
      if(processID==0) say("4sisis>", "Thread ", matejob->id, " completed ", matejob->number," comparisons, joining results");
      int *&result_ref = *matejob->result;
      //					     v--> position + first of the block + first of the thread
      for(size_t i=0; i<matejob->number; i++) result_ref[i + matejob->from_in_result + matejob->from] += result_cp[i];
   pthread_mutex_unlock( matejob->mutex );

   return (void *)0;
}

void nonpareil_count_mates(int *&result, char **&blockA, char **&blockB,
			int fromA, int numberA, int fromB, int numberB,
			int talk, matepar_t matepar){
   // Finally, this is the core of the per-block-per-thread comparisons
   for(int i=0; i<numberA; i++){
      if(processID==0 && talk>0 && i%talk==0) say("4sfs^", "Searching sequences: ", (double)i*100/(double)numberA, "% of the block");
      for(int j=0; j<numberB; j++){
	 if(nonpareil_compare_reads(blockA[i+fromA], blockB[j+fromB], matepar))
	    result[i]++;
      }
   }
   return;
}

bool nonpareil_compare_reads(char *seqA, char *seqB, matepar_t matepar){
   // Vars
   int	lenA = strlen(seqA),
   	lenB = strlen(seqB);

   // Compare
   if(lenA<lenB) return nonpareil_compare_reads_shortfirst(seqA, seqB, lenA, lenB, matepar);
   return nonpareil_compare_reads_shortfirst(seqB, seqA, lenB, lenA, matepar);
}

bool nonpareil_compare_reads_shortfirst(char *seqA, char *seqB, int lenA, int lenB, matepar_t matepar){
   int		min_len = (int)ceil(matepar.overlap * lenA);
   double	disimilarity = 1.0 - matepar.similarity;

   // Compare strands W-W
   for(int i=min_len-lenB; i<=lenA-min_len; i++){
      int	errors = 0,
      		from_seqA = (i<0 ? 0 : i),
      		len_within = (i<0 ? lenB+i : lenA-i),
		max_errors = (int)((disimilarity)*(double)len_within),
		from_seqB = (i<0 ? -1*i : 0);

      for(int j=0; j<len_within; j++){
         if((seqA[j+from_seqA]!=seqB[j+from_seqB] || (matepar.n_as_mismatch && (seqA[j+from_seqA]=='N' || seqB[j+from_seqB]=='N'))) && (++errors>max_errors))
	    goto next_i; // Cannot be 'continue', because it's a nested loop
      }
      return true;

      next_i: ; // To avoid problems with the nested loop.
   }

   // Compare strands C-W
   if(matepar.revcom){
      char	*rc = new char[lenA+1];
      bool	out;

      reverse_complement(rc, seqA);
      matepar.revcom = false;
      out = nonpareil_compare_reads_shortfirst(rc, seqB, lenA, lenB, matepar);
      matepar.revcom = true; // Not really necesary, just in case it's passed as reference or something like that

      free(rc);
      return out;
   }

   return false;
}

void nonpareil_save_mates(int *&result, int no_results, char *file){
   // Vars
   ofstream	fileh;

   fileh.open(file, ios::out);
   if(!fileh.is_open()) error("Cannot open the file", file);

   for(int a=0; a<no_results; a++)
      fileh << result[a] << endl;

   return;
}
