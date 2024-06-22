#!/bin/bash

set -e

for f in fasta fastq ; do
  gzip -cd "test/test.${f}.gz" > "test/test.${f}"
  for t in kmer alignment ; do
    for gz in "" .gz ; do
      echo "=======> TESTING $t $f $gz"
      ./nonpareil -T "$t" -s "test/test.${f}${gz}" -f "$f" -b DELETE -X 50 -v 1
      if [[ $t == alignment && -s nonpareil-mpi ]] ; then
        echo "=======> TESTING MPI $t $f $gz"
        mpirun -np 2 \
          ./nonpareil-mpi -T "$t" -s "test/test.${f}${gz}" -f "$f" -b DELETE \
            -X 50 -v 1
      fi
    done
  done
  rm "test/test.${f}"
done
rm DELETE.np*

