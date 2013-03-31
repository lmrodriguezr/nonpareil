
gpp=g++
enveomics=enveomics/
universal=$(enveomics)universal.o
sqlite=$(enveomics)sqlite.o -lsqlite3
regex=$(enveomics)regex.o
go=$(enveomics)go.o $(regex) $(sqlite)
sequence=$(enveomics)sequence.o
pthread=-lpthread


all:	enveomics nonpareil

enveomics:
	cd $(enveomics) && $(MAKE) all

nonpareil:
	cd $(enveomics) && $(MAKE) nonpareil
	$(gpp) nonpareil.cpp $(universal) $(sequence) $(pthread) $(enveomics)nonpareil_mating.o $(enveomics)nonpareil_sampling.o -o nonpareil

nuc_sampler:
	cd $(enveomics) && $(MAKE) sequence
	$(gpp) nuc_sampler.cpp $(universal) $(sequence) $(pthread) -o utils/nuc_sampler
	
clean:
	cd $(enveomics) && $(MAKE) clean
