
# Makefile for Nonpareil
# @update Dec 16 2013
# @author Luis M. Rodriguez-R <lmrodriguez at gmail dot com>

include globals.mk

enveomics=enveomics/
universal=$(enveomics)universal.o
sqlite=$(enveomics)sqlite.o -lsqlite3
regex=$(enveomics)regex.o
go=$(enveomics)go.o $(regex) $(sqlite)
sequence=$(enveomics)sequence.o
multinode=$(enveomics)multinode.o
pthread=-lpthread
np_objs=$(universal) $(multinode) $(sequence) $(enveomics)nonpareil_mating.o $(enveomics)nonpareil_sampling.o


all:	nonpareil

enveomics:
	cd $(enveomics) && $(MAKE) all

nonpareil-mpi:
	cd $(enveomics) && $(MAKE) $@
	$(mpicpp) nonpareil.cpp $(np_objs) $(pthread) -o $@

nonpareil:
	cd $(enveomics) && $(MAKE) $@
	$(cpp) nonpareil.cpp $(np_objs) $(pthread) -o $@

nuc_sampler:
	cd $(enveomics) && $(MAKE) sequence
	$(cpp) nuc_sampler.cpp $(universal) $(sequence) $(pthread) -o utils/nuc_sampler
	
clean:
	cd $(enveomics) && $(MAKE) clean

install:
	if [ ! -d $(bindir) ] ; then mkdir -p $(bindir) ; fi
	if [ ! -d $(mandir) ] ; then mkdir -p $(mandir) ; fi
	if [ -e nonpareil ] ; then install -m 0755 nonpareil $(bindir)/ ; fi
	if [ -e nonpareil-mpi ] ; then install -m 0755 nonpareil-mpi $(bindir)/ ; fi
	cp docs/_build/man/nonpareil.1 $(mandir)/nonpareil.1
	$(R) CMD INSTALL utils/Nonpareil

