#ifndef SEQREADER_H
#define SEQREADER_H

#include <iostream>
#include <fstream>
#include <unordered_set>
#include <string>
#include <time.h>
#include <stdlib.h>
#include <random>
#include <vector>
#include <string>

using namespace std;

class Sequence {

public:
  Sequence();
  Sequence(string header, string sequence);
  Sequence(string header, string sequence, string qualstr, vector<double> baseProb);
  void setHeader(string header);
  void setSequence(string sequence);
  void setQualStr(string qualstr);
  void setBaseProb(vector<double> baseProb);
  string header, sequence, qualstr;
  vector<double> baseProb;
};

void buildFastqSeq(string header, string sequence, string qual, Sequence &out);
void buildFastaSeq(string header, string sequence);

class SeqReader {

public:
  std::ifstream &ifs;
  string filename;
  std::unordered_set<std::string> randomHeaders;
  bool readNext;
  std::random_device rd;
  std::mt19937_64 gen;
  std::uniform_int_distribution<unsigned long long int> distribution;
  SeqReader(ifstream &ifs);
  void reset();
  virtual size_t readNextSeq(Sequence &out) = 0;
  virtual size_t getRandomSeq(Sequence &out) = 0;
};

class FastqReader: public SeqReader {

public:
  FastqReader(ifstream &ifs);
  size_t readNextSeq(Sequence &out);
  size_t getRandomSeq(Sequence &out);
};

class FastaReader: public SeqReader {
public:
  FastaReader(ifstream &ifs);
  size_t readNextSeq(Sequence &out);
  size_t getRandomSeq(Sequence &out);
};

#endif
