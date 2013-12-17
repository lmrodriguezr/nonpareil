
# Makefile for Nonpareil
# @update Dec 16 2013
# @author Luis M. Rodriguez-R <lmrodriguez at gmail dot com>

include ../globals.mk

mpiflags=-DENVEOMICS_MULTI_NODE
CPPFLAGS=-c -Wall

all: go sequence nonpareil

go:
	$(cpp) universal.cpp $(CPPFLAGS)
	$(cpp) sqlite.cpp $(CPPFLAGS)
	$(cpp) regex.cpp $(CPPFLAGS)
	$(cpp) go.cpp $(CPPFLAGS)

sequence:
	$(cpp) universal.cpp $(CPPFLAGS)
	$(cpp) sequence.cpp $(CPPFLAGS)

nonpareil:
	$(cpp) universal.cpp $(CPPFLAGS)
	$(cpp) multinode.cpp $(CPPFLAGS)
	$(cpp) sequence.cpp $(CPPFLAGS)
	$(cpp) nonpareil_mating.cpp $(CPPFLAGS)
	$(cpp) nonpareil_sampling.cpp $(CPPFLAGS)

nonpareil-mpi:
	$(mpicpp) $(mpiflags) universal.cpp $(CPPFLAGS)
	$(mpicpp) $(mpiflags) multinode.cpp $(CPPFLAGS)
	$(mpicpp) $(mpiflags) sequence.cpp $(CPPFLAGS)
	$(mpicpp) $(mpiflags) nonpareil_mating.cpp $(CPPFLAGS)
	$(mpicpp) $(mpiflags) nonpareil_sampling.cpp $(CPPFLAGS)
	
clean:
	rm *.o

