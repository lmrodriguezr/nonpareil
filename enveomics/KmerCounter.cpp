#include "Hash.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include "KmerCounter.h"
#include "universal.h"

using namespace std;

KmerCounter::KmerCounter(References &references, FastqReader &metaReader, string cntFile) {
  counting(references,metaReader);
  prepCounts(references);
  saveCounts(cntFile);
}

KmerCounter::KmerCounter(References &references, FastaReader &metaReader, string cntFile) {
  counting(references,metaReader);
  prepCounts(references);
  saveCounts(cntFile);
}

void KmerCounter::counting(References &references, FastqReader &metaReader) {
  Sequence temp;
  Hash hasher(references.ksize);
  unsigned long long int hashcode;
  while(metaReader.readNextSeq(temp) != (size_t)(-1)) {
    if(temp.sequence.length() < references.ksize)
        error("Reads are required to have a minimum length of kmer size");
    hasher.intialize(temp.sequence);
    this->totalSeqs++;
    this->totalLength = this->totalLength + temp.sequence.length();
    hashcode = 0;
    while(hasher.nextHash(hashcode) != (size_t)(-1)) {
      references.update(hashcode);
    }
  }
}

void KmerCounter::counting(References &references, FastaReader &metaReader) {
  Sequence temp;
  Hash hasher(references.ksize);
  unsigned long long int hashcode;
  while(metaReader.readNextSeq(temp) != (size_t)(-1)) {
    cout << temp.sequence << temp.sequence.length() << endl;
    if(temp.sequence.length() < references.ksize)
        error("Reads are required to have a minimum length of kmer size");
    hasher.intialize(temp.sequence);
    this->totalSeqs++;
    this->totalLength = this->totalLength + temp.sequence.length();
    hashcode = 0;
    while(hasher.nextHash(hashcode) != (size_t)(-1)) {
      references.update(hashcode);
    }
  }
}

void KmerCounter::prepCounts(References &references) {
  int count = 0;
  for(size_t i = 0; i < references.refKmers.size(); i++) {
    count = references.refKmerMap[references.refKmers[i]] + references.refKmerMap[references.refRevComKmers[i]];
    this->countTable.push_back(count);
  }

  size_t j = 0;
  for(size_t i = 0; i < nearbyint(references.totalErrKmers); i++) {
    while(true){
      if(countTable.size() <= j)
          break;
      if(countTable[j] == 1) {
        countTable.erase(countTable.begin()+j);
        break;
      }
      j++;
    }
  }
}

void KmerCounter::saveCounts(string file) {
  ofstream	fileh;

  fileh.open(file, ios::out);
  if(!fileh.is_open()) error("Cannot open the file", file);

  for(size_t a=0; a<this->countTable.size(); a++)
     fileh << countTable[a] << endl;
}

unsigned int KmerCounter::getTotalSeqs() {
  return totalSeqs;
}

int KmerCounter::getTotalQSeqs() {
  return countTable.size();
}

double KmerCounter::getAvgLen() {
  double temp;
  temp = double(this->totalLength)/double(this->totalSeqs);
  return temp;
}

void KmerCounter::getCounts(int *&result) {
  //result = new int[this->countTable.size()];

  for(size_t i = 0; i < countTable.size(); i++) {
    result[i] = countTable[i];
  }
}
