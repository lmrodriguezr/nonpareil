Redundancy
==========

If you have large files (>1Gb) and access to a cluster, take a look at :doc:`mpi`.

For the impatient
-----------------

Even if you're in a hurry, taking a look at :doc:`preprocess` is very important. If you already did, you can simply run::

    nonpareil -s reads.fa -b output

Where ``reads.fa`` is the file containing the trimmed single reads, and ``output`` is the prefix
of the output files to be created.

Mandatory options
-----------------
   -s <str>   Path to the (input) file containing the sequences.  This is lowercase S.

Common options
--------------
   -f <str>   The format of the sequences.  Can be 'fasta' or 'fastq'.  By default: 'fasta'
   -b <str>   Path to the prefix for all the output files.  Replaces the options: -a, -C, -l, and -o; generating files
              with the suffixes .npa, npc, .npl, and .npo, respectively, unless explicitly set.
   -i <num>   Interval between sampling portions.  By default: 0.01.
   -n <int>   Number of sub-samples to generate per point.  If it is not a multiple of the number of threads (see -t),
              it is rounded to the next (upper) multiple.  By default: 1024.
   -L <num>   Minimum overlapping percentage of the aligned region on the largest sequence. The similarity (see -S) is
              evaluated for the aligned region only.  By default: 50.
   -X <int>   Maximum number of reads to use as query.  This is capital X.  By default, 1,000 reads.
   -R <int>   Maximum RAM usage in Mib.  Ideally this value should be larger than the sequences to analyze (discarding
              non-sequence elements like headers or quality).  This is particularly important when running in multiple
              cores (see -t).  This value is approximated.  By default 1024.
              Maximum value in this version: 4194303
   -t <int>   Number of threads.  Highest efficiency when the number of sub-samples (see -n) is multiple of the number
              of threads.  By default: 2.
   -v <int>   Verbosity level, for debugging purposes.  By default 7.  This is lowercase V.
   -V         Show version information and exit.  This is uppercase V.
   -h         Display this message and exit.

Additional options
------------------
**Input/Output**
   -a <str>   Path to the (output) file where all data must be saved.  This report is not created by default.  See the
              OUTPUT section.
   -C <str>   Path to the (output) file where the mating vector is to be saved.  This is a capital C.
   -F         Report the sampled portions as a fraction of the library instead of the number of reads.  See -a, -o and
              the OUTPUT section.
   -l <str>   Path to the (output) file where the log of the run must be saved. By default the log is sent only to the
              STDERR.  If set, the log is sent to both the STDERR and the log file.
   -o <str>   Path to the (output) file where summary is to be saved.   By default the summary is sent to stdout (same
              behavior as using a dash '-').  If an empty string '' is provided, does not produce the summary. See the
              OUTPUT section.
   
**Sampling**
   -m <num>   Minimum value of sampling portion.  By default: 0.
   -M <num>   Maximum value of sampling portion.  By default: 1.
   -d <num>   Take this fraction of the total library every sampling point (logarithmic sampling, not linear).  If set
              to zero, logarithmic sub-sampling is disabled (default). Recommended value: 0.7. EXPERIMENTAL CODE.
   
**Mating**
   -c         Do not use reverse-complement.  This is useful for single stranded sequences data (like RNA).  This is a
              lowercase C.
   -N         Treat Ns as mismatches.  By default, Ns (unknown nucleotides) match any nucleotide (even another N).
   -S <num>   Similarity threshold to group two reads together.   Reducing this option will increase sensitivity while
              increasing running time.  This is uppercase S.
   -x <num>   Probability of taking a sequence into account as query for the construction of the curve.  Higher values
              reduce accuracy but increase speed.  This is lower case x.  If set, overides -X.

**Misc**
   -A         Autoadjust parameters and re-run.  Evaluates the results looking for common problems, adjusts parameters
              and re-run the analyses.  THIS IS EXPERIMENTAL CODE.
   -r <int>   Random generator seed.  By default current time.

Input
-----
Sequences must be in FastA or FastQ format. See :doc:`preprocess`.

Output
------
Redundancy summary: ``.npo`` file
   Tab-delimited file with six columns. The first column indicates the sequencing effort (in number of reads), and the
   remaining columns indicate the summary of the distribution of redundancy (from the replicates, 1,024 by default) at
   the given sequencing effort. These five columns are: average redundancy, standard deviation, quartile 1, median
   (quartile 2), and quartile 3.

Redundancy values: ``.npa`` file
   Tab-delimited file with three columns. Similar to the .npo files, it contains information about the redundancy at
   each sequencing effort, but it provides ALL the results from the replicates, not only the summary at each point. The
   first column indicates the sequencing effort (as a fraction of the dataset), the second column indicates the ID of
   the replicate (a number used only to introduce some controlled noise in plots), and the third column indicates the
   estimated redundancy value.

Mates distribution: ``.npc`` file
   Raw list with the number of reads in the dataset matching a query read. A set of query reads is randomly drawn by
   Nonpareil (1,000 by default), and compared against all reads in the dataset. Each line on this file corresponds to a
   query read (the order is not important). We have seen certain correspondance between these numbers and the distribution
   of abundances in the community (compared, for example, as rank-abundance plots), but this file is provided only for
   quality-control purposes and comparisons with other tools.

Log: ``.npl`` file
   A verbose log of internal Nonpareil processing. The number to the left (inside squared brackets) indicate the CPU time
   (in minutes). This file also provide quality assessment of the Nonpareil run (automated consistency evaluation). Ideally,
   the last line should read "Everything seems correct". Otherwise, it suggests alternative parameters that may improve the
   estimation.

