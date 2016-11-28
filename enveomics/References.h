#ifndef REFERENCES_H
#define REFERENCES_H

#include <unordered_map>
#include <vector>
#include <string>
#include "SeqReader.h"
using namespace std;

class References {

public:
  double totalErrKmers;
  size_t refSize, ksize;
  vector<unsigned long long int> refKmers;
  vector<unsigned long long int> refRevComKmers;
  unordered_map<unsigned long long int, int> refKmerMap;
  References(FastqReader &fastqReader, int refSize, int ksize);
  References(FastaReader &fastaReader, int refSize, int ksize);
  References(FastaReader &fastaReader, int ksize, bool alt_query);
  void update(unsigned long long int hashcode);

private:
  void intializeReferences(FastqReader &fastqReader);
  void intializeReferences(FastaReader &fastaReader);
  void intializeReferences(FastaReader &fastaReader, bool alt_query);
};

#endif
