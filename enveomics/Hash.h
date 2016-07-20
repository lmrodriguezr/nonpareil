#ifndef HASH_H
#define HASH_H

#include <string>

using namespace std;


size_t getHashCode(string s, unsigned long long int &hashcode);

class Hash {
public:
  Hash(int k);
  void intialize(string seq);
  size_t nextHash(unsigned long long int &hashcode);

private:
  int k;
  string seq;
  int len;
  int pos;
  int size;
  unsigned long long int kmask;
};

#endif
