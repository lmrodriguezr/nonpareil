
gpp=g++

all:
	$(gpp) universal.cpp -c
	$(gpp) sqlite.cpp -c
	$(gpp) regex.cpp -c
	$(gpp) go.cpp -c

go:
	$(gpp) universal.cpp -c
	$(gpp) sqlite.cpp -c
	$(gpp) regex.cpp -c
	$(gpp) go.cpp -c

sequence:
	$(gpp) universal.cpp -c
	$(gpp) sequence.cpp -c

nonpareil:
	$(gpp) universal.cpp -c
	$(gpp) sequence.cpp -c
	$(gpp) nonpareil_mating.cpp -c
	$(gpp) nonpareil_sampling.cpp -c

clean:
	rm *.o

