#ifndef KMER_COUNTER_H
#define KMER_COUNTER_H

#include "References.h"

using namespace std;

class KmerCounter {

public:
  KmerCounter(References &references, FastqReader &metaReader, string cntFile);
  KmerCounter(References &references, FastaReader &metaReader, string cntFile);
  int getCountTable(int *&result);
  unsigned int getTotalSeqs();
  double getAvgLen();
  void getCounts(int *&result);
  int getTotalQSeqs();
private:
  vector<int> countTable;
  unsigned long long int totalSeqs;
  unsigned long long int totalLength;
  void counting(References &references, FastqReader &metaReader);
  void counting(References &references, FastaReader &metaReader);
  void prepCounts(References &references);
  void saveCounts(string file);
};

#endif
