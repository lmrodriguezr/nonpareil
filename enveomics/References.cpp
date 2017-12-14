#include "References.h"
#include "Hash.h"
#include "sequence.h"
#include "universal.h"

using namespace std;

References::References(FastqReader &fastqReader, int refSize, int ksize) {
  this->refSize = refSize;
  this->ksize = ksize;
  this->totalErrKmers = 0;
  intializeReferences(fastqReader);
}

References::References(FastaReader &fastaReader, int refSize, int ksize) {
  this->refSize = refSize;
  this->ksize = ksize;
  //this->totalErrKmers = 0.014 * refSize;
  this->totalErrKmers = 0;
  intializeReferences(fastaReader);
}

References::References(FastaReader &fastaReader, int ksize, bool alt_query) {
  this->ksize = ksize;
  this->refSize = 0;
  this->totalErrKmers = 0;
  intializeReferences(fastaReader, alt_query);
}

void References::intializeReferences(FastaReader &fastaReader, bool alt_query) {
  Sequence temp;
  string kmer, revkmer;
  int flag;
  unsigned long long int hashcode;

  while(fastaReader.readNextSeq(temp) != (size_t)(-1)) {
    this->refSize++;
    if(temp.sequence.length() < this->ksize)
        error("Reads are required to have a minimum length of kmer size");
    kmer = temp.sequence.substr(0,this->ksize);
    flag = getHashCode(kmer, hashcode);
    if(flag == -1) {
      this->refSize--;
      continue;
    }
    refKmerMap[hashcode] = 0;
    refKmers.push_back(hashcode);

    reverse_complement(revkmer, kmer);
    getHashCode(revkmer, hashcode);
    refKmerMap[hashcode] = 0;
    refRevComKmers.push_back(hashcode);
  }
  //this->totalErrKmers = 0.014 * refSize;
}

void References::intializeReferences(FastqReader &fastqReader) {
  Sequence temp;
  string kmer, revkmer;
  int flag;
  double kerr;
  unsigned long long int hashcode;
  size_t i = 0;
  // deal with shorter reads then kmer length
  for(i=0;i<this->refSize;i++) {
    fastqReader.getRandomSeq(temp);
    //Hashcode for forward kmer
    if(temp.sequence.length() < this->ksize)
      error("Reads are required to have a minimum length of kmer size");
    kmer = temp.sequence.substr(0,this->ksize);
    flag = getHashCode(kmer,hashcode);
    if(flag == -1) {
        i--;
        continue;
    }
    refKmerMap[hashcode] = 0;
    refKmers.push_back(hashcode);

    //Hashcode reverse complement kmer
    reverse_complement(revkmer, kmer);
    getHashCode(revkmer,hashcode);
    refKmerMap[hashcode] = 0;
    refRevComKmers.push_back(hashcode);
    kerr = 1.0;
    for(size_t j = 0; j < ksize; j++) {
      kerr = kerr * (1.0 - temp.baseProb[j]);
    }
    this->totalErrKmers = this->totalErrKmers + (1-kerr);
  }
}

void References::intializeReferences(FastaReader &fastaReader) {
  Sequence temp;
  string kmer, revkmer;
  int flag;
  unsigned long long int hashcode;
  size_t i = 0;
  for(i=0;i<this->refSize;i++) {
    fastaReader.getRandomSeq(temp);
    if(temp.sequence.length() < this->ksize)
        error("Reads are required to have a minimum length of kmer size");
    kmer = temp.sequence.substr(0,this->ksize);
    flag = getHashCode(kmer, hashcode);
    if(flag == -1) {
      i--;
      continue;
    }
    refKmerMap[hashcode] = 0;
    refKmers.push_back(hashcode);

    reverse_complement(revkmer, kmer);
    getHashCode(revkmer, hashcode);
    refKmerMap[hashcode] = 0;
    refRevComKmers.push_back(hashcode);
  }
}

void References::update(unsigned long long int hashcode){
  if(this->refKmerMap.count(hashcode) > 0) {
    refKmerMap[hashcode]++;
  }
}
