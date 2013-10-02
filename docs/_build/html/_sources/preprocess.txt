Preprocess the reads
====================

Nonpareil expects that the sequencing error is always well below 5%, so we suggest using an expected error cutoff of 1%
(i.e., Q>20, or 1 error in 100 nucleotides). We recommend to perform this task using SolexaQA.

Ideally, the reads should be in FastA format (althought Nonpareil can read FastQ). To transform FastQ into FastA, you
can simply use::

    # Input: reads.fastq
    # Output: reads.fasta
    cat reads.fastq | paste - - - - | awk 'BEGIN{FS="\t"}{print ">"substr($1,2)"\n"$2}' > reads.fasta


Also, if you have paired-end reads, you should use only one sister read per pair in Nonpareil. If you have them interposed
in the same file, you can separate them using FastA.split.pl_::

    # Input: reads.fasta
    # Output: reads.1.fa and reads.2.fa
    FastA.split.pl reads.fasta reads 2

.. _FastA.split.pl: https://github.com/lmrodriguezr/enveomics/blob/master/Scripts/FastA.split.pl

