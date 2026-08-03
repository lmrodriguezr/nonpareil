Redundancy
==========

If you have large files (>1Gb) and access to a cluster, take a look at
:doc:`mpi`.

For the impatient
-----------------

Even if you're in a hurry, taking a look at :doc:`preprocess` is very important.
If you already did, you can simply run::

    # fastq is recommended for kmer algorithm (faster, short reads)
    nonpareil -s reads.fa -T kmer -f fastq -b output
    nonpareil -s reads.fa -T kmer -f fasta -b output

    # fasta is recommended for alignment algorithm (slower, short reads)
    nonpareil -s reads.fa -T alignment -f fasta -b output
    nonpareil -s reads.fa -T alignment -f fastq -b output

    # gzipped files are also supported (with .gz extension)
    nonpareil -s reads.fa.gz -T kmer -f fastq -b output

    # Usearch-based nonpareil with FastA input (long reads)
    nonpareil -s reads.fa -T usearch -f fasta -b output

    # general information about Nonpareil
    nonpareil -h # A help message
    nonpareil -V # Current version of Nonpareil
    nonpareil -B # Bibliography, how to cite Nonpareil


The input files in the examples above are ``reads.fa``, ``reads.fq``, or their
gzipped versions, and ``output`` is the prefix of the output files to be
created.

Mandatory options
-----------------
   -s <str>   Path to the (input) file containing the sequences.
              Gzipped files are supported with .gz extension
   -T <str>   nonpareil algorithm. can be 'kmer' or 'alignment'
   -f <str>   The format of the sequence: 'fasta' or 'fastq'

Common options
--------------
   -b <str>   Path to the prefix for all the output files.  Replaces the
              options: -a, -C, -l, and -o; generating files with the suffixes
              .npa, npc, .npl, and .npo, respectively, unless explicitly set
   -X <int>   Maximum number of reads to use as query.  This is capital X.  By
              default, 1,000 reads
   -k <int>   kmer size. You can increase kmer size to increase sensitivity.  By
              default: 24.  Only applies to ``-T kmer``
   -n <int>   Number of sub-samples to generate per point.  If it is not a
              multiple of the number of threads (see -t), it is rounded to the
              next (upper) multiple.  By default: 1024
   -L <num>   Minimum overlapping percentage of the aligned region on the
              largest sequence. The similarity (see -S) is evaluated for the
              aligned region only.  By default: 50
   -p <int>   Minimum read length (base pairs).  By default: 50 for ``-T kmer``
              and ``-T alignment``, 1000 for ``-T usearch``

   -R <int>   Maximum RAM usage in Mib.  Ideally this value should be larger
              than the sequences to analyze (discarding non-sequence elements
              like headers or quality).  This is particularly important when
              running in multiple cores (see -t).  This value is approximated.
              By default 1024. Maximum value in this version: 4194303
   -t <int>   Number of threads.  Highest efficiency when the number of
              sub-samples (see -n) is multiple of the number of threads.  By
              default: 2
   -v <int>   Verbosity level, for debugging purposes.  By default 7.  This is
              lowercase V
   -r <int>   Random generator seed.  By default current time
   -V         Show version information and exit.  This is uppercase V
   -B         Show bibliography (citation) entries and exit
   -h         Display this message and exit

Additional options
------------------
**Input/Output**
   -a <str>   Path to the (output) file where all data must be saved.  This
              report is not created by default.  See the OUTPUT section
   -C <str>   Path to the (output) file where the mating vector is to be saved.
              This is a capital C
   -F         Report the sampled portions as a fraction of the library instead
              of the number of reads.  See -a, -o and the OUTPUT section
   -l <str>   Path to the (output) file where the log of the run must be saved.
              By default the log is sent only to the STDERR.  If set, the log is
              sent to both the STDERR and the log file
   -o <str>   Path to the (output) file where summary is to be saved.  By
              default the summary is sent to stdout (same behavior as using a
              dash '-').  If an empty string '' is provided, does not produce
              the summary. See the OUTPUT section
   -q <str>   Path to the (input) file containing a second dataset to be used as
              query, for dataset comparisons.  This option is currently
              experimental
   -K         Keep the sequence index (the ``.enve-seq``/``.enve-nam`` files
              generated next to the input) after the run finishes, instead of
              removing it. These files can always be reused across the MPI
              processes of a single run regardless of this option; ``-K`` only
              affects whether they also survive to be reused by a *separate*,
              later run against the same input (skipping re-indexing). Has no
              effect if the input's directory isn't writable, since in that
              case the index is built in a temporary directory that is always
              removed at the end of the run. This is uppercase K

**Sampling**
   -d <num>   Subsample iteratively applying this factor to the number of reads,
              resulting in logarithmic subsampling. Use -d 0 to fall back to
              linear sampling, controlled by -m, -M, & -i (this was the default
              before v2.4). By default: 0.7
   -m <num>   Minimum value of sampling portion.  By default: 0
   -M <num>   Maximum value of sampling portion.  By default: 1
   -i <num>   Interval between sampling portions. By default: 0.01

**Mating**
   -c         Do not use reverse-complement.  This is useful for single stranded
              sequences data (like RNA). Only applies to ``-T alignment``.  This
              is a lowercase C
   -H <num>   Hash size. The number of slots to use in the read indexing step.
              Ideally, this should be 1.5-2 times the total number of reads. By
              default: 0, meaning that slots are set by usearch. Only applies to
              ``-T usearch``
   -N         Treat Ns as mismatches.  By default, Ns (unknown nucleotides)
              match any nucleotide (even another N).  Only applies to
              ``-T alignment``
   -S <num>   Similarity threshold to group two reads together.  Reducing this
              option will increase sensitivity while increasing running time.
              Only applies to ``-T alignment``. This is uppercase S
   -x <num>   Probability of taking a sequence into account as query for the
              construction of the curve.  Higher values reduce accuracy but
              increase speed.  This is lower case x.  If set, overides -X

**Misc**
   -A         Autoadjust parameters and re-run.  Evaluates the results looking
              for common problems, adjusts parameters and re-run the analyses.
              THIS IS EXPERIMENTAL CODE

Input
-----
Sequences must be in FastA or FastQ format. The input file can be gzipped
provided it has ``.gz`` extension.  See :doc:`preprocess`.

Output
------
Redundancy summary: ``.npo`` file
   Tab-delimited file with six columns. The first column indicates the
   sequencing effort (in number of reads), and the remaining columns indicate
   the summary of the distribution of redundancy (from the replicates, 1,024 by
   default) at the given sequencing effort. These five columns are: average
   redundancy, standard deviation, quartile 1, median (quartile 2), and
   quartile 3.

Redundancy values: ``.npa`` file
   Tab-delimited file with three columns. Similar to the .npo files, it contains
   information about the redundancy at each sequencing effort, but it provides
   ALL the results from the replicates, not only the summary at each point. The
   first column indicates the sequencing effort (as a fraction of the dataset),
   the second column indicates the ID of the replicate (a number used only to
   introduce some controlled noise in plots), and the third column indicates the
   estimated redundancy value.

Mates distribution: ``.npc`` file
   Raw list with the number of reads in the dataset matching a query read. A set
   of query reads is randomly drawn by Nonpareil (1,000 by default), and
   compared against all reads in the dataset. Each line on this file corresponds
   to a query read (the order is not important). We have seen certain
   correspondance between these numbers and the distribution of abundances in
   the community (compared, for example, as rank-abundance plots), but this file
   is provided only for quality-control purposes and comparisons with other
   tools.

Log: ``.npl`` file
   A verbose log of internal Nonpareil processing. The number to the left
   (inside squared brackets) indicates the wall-clock time elapsed since the
   run started (in minutes). This file also provides quality assessment of
   the Nonpareil run (automated consistency evaluation). Ideally, the last
   line should read "Everything seems correct".
   Otherwise, it suggests alternative parameters that may improve the
   estimation.

