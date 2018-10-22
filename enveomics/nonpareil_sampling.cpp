// nonpareil_sampling - Part of the nonpareil package
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
#include <vector>
#include <algorithm>

#include "universal.h"
// #include "multinode.h"
#include "sequence.h"
#include "nonpareil_sampling.h"

//extern int processID;
//extern int processes;

#define LARGEST_LABEL	128
#define LARGEST_LINE	2048

using namespace std;

int nonpareil_sample_portion(double *&result, int threads, samplepar_t samplepar){
   // Vars
   if(samplepar.replicates<threads) threads = samplepar.replicates;
   samplejob_t		samplejob[threads];
   pthread_t		thread[threads];
   pthread_mutex_t	mutex = PTHREAD_MUTEX_INITIALIZER;
   int			rc, samples_per_thr, launched_replicates=0;

   // Blank result
   result = new double[samplepar.replicates];
   for(int a=0; a<samplepar.replicates; a++) result[a] = 0.0;

   // Set sample
   //if(processID==0)
   say("4sfs^", "Sampling at ", samplepar.portion*100, "%");
   samples_per_thr = (int)ceil((double)samplepar.replicates/(double)threads);
   threads = (int)ceil((double)samplepar.replicates/samples_per_thr);

   // Launch samplings
   for(int thr=0; thr<threads; thr++){
      samplejob[thr].id = thr;
      samplejob[thr].from_in_result = thr*samples_per_thr; // Zero-based index to start at in the results array
      samplejob[thr].number = samplejob[thr].from_in_result+samples_per_thr > samplepar.replicates ? samplepar.replicates-samplejob[thr].from_in_result : samples_per_thr;
      samplejob[thr].samplepar = samplepar;
      samplejob[thr].result = &result;
      samplejob[thr].mutex = &mutex;
      launched_replicates += samplejob[thr].number;

      if((rc=pthread_create(&thread[thr], NULL, &nonpareil_sample_portion_thr, (void *)&samplejob[thr]  )))
	 error("Thread creation failed", (char)rc);
   }

   // Gather jobs
   for(int thr=0; thr<threads; thr++){
      //if(thr%processes == processID)
      pthread_join(thread[thr], NULL);
   }

   // Return
   return launched_replicates;
}

void *nonpareil_sample_portion_thr(void *samplejob_ref){
   // Vars
   samplejob_t	*samplejob = (samplejob_t *)samplejob_ref;
   int		*&mates_ref = *samplejob->samplepar.mates, n;
   double	*result_cp = new double[samplejob->number], p, p_gt_0;

   // Sample
   for(int i=0; i<samplejob->number; i++){
      int	found=0, sample_size=0;
      // For each read with known number of mates
      for(int j=0; j<samplejob->samplepar.mates_size; j++){
         // Include the read in the sample?
	 if(((double)rand()/(double)RAND_MAX) < samplejob->samplepar.portion){
	    sample_size++;
	    // Does the sample contain at least one mate?
	    // if((double)rand()/(double)RAND_MAX < 1.0-pow(1-samplejob->samplepar.portion, mates_ref[j]-1)) found++;
	    n=(int)floor((samplejob->samplepar.total_reads)*samplejob->samplepar.portion)-1; // -1 because we exclude the query read from the sample
	    if(n<0) n=0;
	    p=(double)(mates_ref[j]-1.0)/(samplejob->samplepar.total_reads-1.0); // Again, -1 in both terms because of the query read
	    p_gt_0 = 1.0 - pow(1.0-p, n);
	    //cout << "n=" << n << ", p=" << p << "P(X>0)=" << p_gt_0 << endl;
	    if((double)rand()/(double)RAND_MAX < p_gt_0) found++;
	 }
      }
      result_cp[i] = sample_size==0 ? 0.0 : (double)found/(double)sample_size;
   }

   // Transfer results to the external vector
   pthread_mutex_lock( samplejob->mutex );
      double *&result_ref = *samplejob->result;
      //					        v--> position + first of the thread
      for(int i=0; i<samplejob->number; i++) result_ref[i + samplejob->from_in_result] = result_cp[i];
   pthread_mutex_unlock( samplejob->mutex );

   return (void *)0;
}

sample_t nonpareil_sample_summary(double *&sample_result, int sample_number, char *alldata, char *outfile, samplepar_t samplepar){
   // Vars
   double		x2=0;
   vector<double>	dataPoints;
   FILE			*alldatah, *summaryh;
   bool			reportAllData=false;
   int			reportSummary=0;
   char			*text, *label, *sep, *header=(char*)"";
   sample_t		s;

   s.portion = samplepar.portion;

   // File handlers
   if(alldata && strlen(alldata)>0){
      alldatah = fopen(alldata,
        (samplepar.portion==samplepar.portion_min ? "w" : "a+"));
      if(alldatah==NULL) error("I cannot write on all-data file", alldata);
      reportAllData=true;
   }
   if(strlen(outfile)>0 && strcmp(outfile, "-")!=0){
      summaryh = fopen(outfile,
        (samplepar.portion==samplepar.portion_min ? "w" : "a+"));
      if(summaryh==NULL) error("I cannot write on summary file", outfile);
      reportSummary=1;
   }else if(strlen(outfile)>0){
      reportSummary=2;
   }

   if(samplepar.portion==samplepar.portion_min){
      header = new char[LARGEST_LINE];
      if(samplepar.type == 1) { // Kernel: Alignment
        sprintf(header,
          "# @impl: Nonpareil\n# @version: %.2f\n# @maxL: %d\n# @L: %.3f\n# @R: %llu\n# @overlap: %.2f\n# @divide: %.2f\n",
          samplepar.np_version,		// @version
          samplepar.max_read_len,	// @maxL
          samplepar.avg_read_len,	// @L
          samplepar.total_reads,	// @R
          samplepar.seq_overlap*100.0,	// @overlap
          samplepar.divide);		// @divide
      }else if(samplepar.type == 2) { // Kernel: Kmer
        sprintf(header,
          "# @impl: Nonpareil\n# @ksize: %d\n# @version: %.2f\n# @L: %.3f\n# @AL: %.3f\n# @R: %llu\n# @overlap: %.2f\n# @divide: %.2f\n",
          samplepar.k,			// @ksize
          samplepar.np_version,		// @version
          samplepar.avg_read_len,	// @L
          samplepar.adj_avg_read_len,	// @AL
          samplepar.total_reads,	// @R
          samplepar.seq_overlap*100.0,	// @overlap
          samplepar.divide);		// @divide
      }
   }

   if(sample_number>0){
      // Average & SD
      s.avg = s.sd = x2 = 0.0;
      for(int a=0; a<sample_number; a++){
	 s.avg += (double)sample_result[a];
	 x2  += pow((double)sample_result[a], 2.0);
	 dataPoints.push_back(sample_result[a]);
	 if(reportAllData)
	   fprintf(alldatah, "%.2f\t%d\t%.10f\n",
	     samplepar.portion, a, sample_result[a]);
      }
      fclose(alldatah);
      s.avg /= (double)sample_number;
      x2    /= (double)sample_number;
      s.sd   = sqrt(x2-pow(s.avg, 2.0));

      // Quartiles
      vector<double>::iterator	firstDatum = dataPoints.begin();
      vector<double>::iterator	lastDatum = dataPoints.end();
      vector<double>::iterator	q1Datum = firstDatum + (lastDatum - firstDatum) / 4;
      vector<double>::iterator	q2Datum = firstDatum + (lastDatum - firstDatum) / 2;
      vector<double>::iterator	q3Datum = firstDatum + (lastDatum - firstDatum) * 3 / 4;
      nth_element(firstDatum, q1Datum, lastDatum);
      nth_element(firstDatum, q2Datum, lastDatum);
      nth_element(firstDatum, q3Datum, lastDatum);
      s.q1 = *q1Datum;
      s.q2 = *q2Datum;
      s.q3 = *q3Datum;
   }else{
      s.avg = 0;
      s.sd  = 0;
      s.q1  = 0;
      s.q2  = 0;
      s.q3  = 0;
   }

   // Report summary
   if(reportSummary==0) return s;
   label = new char[LARGEST_LABEL];
   text = new char[LARGEST_LINE];
   sep = (char *)"\t";
   sprintf(label, (samplepar.portion_as_label ? "%.6f" : "%.0f"), samplepar.portion*(samplepar.portion_as_label? 1.0 : samplepar.total_reads ));
   sprintf(text, "%s%s%.5f%s%.5f%s%.5f%s%.5f%s%.5f",
   		label, sep, s.avg, sep, s.sd, sep, s.q1, sep, s.q2, sep, s.q3);
   if(reportSummary==1) {
      fprintf(summaryh, "%s%s\n", header, text);
      fclose(summaryh);
   } else if(reportSummary==2) printf("%s%s\n", header, text);

   return s;
}
