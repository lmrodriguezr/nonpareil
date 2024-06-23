
# Makefile for Nonpareil
# @update Dec 16 2013
# @author Luis M. Rodriguez-R <lmrodriguez at gmail dot com>

include ../globals.mk

all: go sequence nonpareil

go:
	$(cpp) universal.cpp -c $(CPPFLAGS)
	$(cpp) sqlite.cpp -c $(CPPFLAGS)
	$(cpp) regex.cpp -c $(CPPFLAGS)
	$(cpp) go.cpp -c $(CPPFLAGS)

sequence:
	$(cpp) universal.cpp -c $(CPPFLAGS)
	$(cpp) sequence.cpp -c $(CPPFLAGS)

nonpareil:
	$(cpp) universal.cpp -c $(CPPFLAGS)
	$(cpp) multinode.cpp -c $(CPPFLAGS)
	$(cpp) sequence.cpp -c $(CPPFLAGS)
	$(cpp) nonpareil_mating.cpp -c $(CPPFLAGS)
	$(cpp) nonpareil_sampling.cpp -c $(CPPFLAGS)
	$(cpp) SeqReader.cpp -c $(CPPFLAGS)
	$(cpp) References.cpp -c $(CPPFLAGS)
	$(cpp) Hash.cpp -c $(CPPFLAGS)
	$(cpp) KmerCounter.cpp -c $(CPPFLAGS)

nonpareil-mpi:
	$(mpicpp) $(mpiflags) universal.cpp -c $(CPPFLAGS)
	$(mpicpp) $(mpiflags) multinode.cpp -c $(CPPFLAGS)
	$(mpicpp) $(mpiflags) sequence.cpp -c $(CPPFLAGS)
	$(mpicpp) $(mpiflags) nonpareil_mating.cpp -c $(CPPFLAGS)
	$(mpicpp) $(mpiflags) nonpareil_sampling.cpp -c $(CPPFLAGS)
	$(mpicpp) $(mpiflags) SeqReader.cpp -c $(CPPFLAGS)
	$(mpicpp) $(mpiflags) References.cpp -c $(CPPFLAGS)
	$(mpicpp) $(mpiflags) Hash.cpp -c $(CPPFLAGS)
	$(mpicpp) $(mpiflags) KmerCounter.cpp -c $(CPPFLAGS)

clean:
	rm -f *.o
