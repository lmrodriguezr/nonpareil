#include "Hash.h"

unsigned long long int lookupTable[128][2];

void setupLookup() {
  lookupTable[(unsigned char)'A'][0] = 0;
  lookupTable[(unsigned char)'A'][1] = -1;
  lookupTable[(unsigned char)'a'][0] = 0;
  lookupTable[(unsigned char)'a'][1] = -1;
  lookupTable[(unsigned char)'G'][0] = 1;
  lookupTable[(unsigned char)'G'][1] = -1;
  lookupTable[(unsigned char)'g'][0] = 1;
  lookupTable[(unsigned char)'g'][1] = -1;
  lookupTable[(unsigned char)'C'][0] = 2;
  lookupTable[(unsigned char)'C'][1] = -1;
  lookupTable[(unsigned char)'c'][0] = 2;
  lookupTable[(unsigned char)'c'][1] = -1;
  lookupTable[(unsigned char)'T'][0] = 3;
  lookupTable[(unsigned char)'T'][1] = -1;
  lookupTable[(unsigned char)'t'][0] = 3;
  lookupTable[(unsigned char)'t'][1] = -1;
}

size_t getHashCode(string s, unsigned long long int &hashcode) {
  if(lookupTable[(unsigned char)'G'][0] != 1)
    setupLookup();

  if(s.length() > 32) {
    return -1;
  }

  int flag = 1;
  hashcode = 0;

  for(int i = 0; i < s.length(); i++) {
    hashcode = hashcode << 2;
    hashcode |= lookupTable[(unsigned char)s[i]][0];
    hashcode &= lookupTable[(unsigned char)s[i]][1];
    flag &= lookupTable[(unsigned char)s[i]][1];
    if(flag == 0)
      return -1;
  }
  return 0;
}

Hash::Hash(int k) {
  this->k = k;
  this->pos = 0;
  this->kmask = 0;
  lookupTable[(unsigned char)'A'][0] = 0;
  lookupTable[(unsigned char)'A'][1] = -1;
  lookupTable[(unsigned char)'a'][0] = 0;
  lookupTable[(unsigned char)'a'][1] = -1;
  lookupTable[(unsigned char)'G'][0] = 1;
  lookupTable[(unsigned char)'G'][1] = -1;
  lookupTable[(unsigned char)'g'][0] = 1;
  lookupTable[(unsigned char)'g'][1] = -1;
  lookupTable[(unsigned char)'C'][0] = 2;
  lookupTable[(unsigned char)'C'][1] = -1;
  lookupTable[(unsigned char)'c'][0] = 2;
  lookupTable[(unsigned char)'c'][1] = -1;
  lookupTable[(unsigned char)'T'][0] = 3;
  lookupTable[(unsigned char)'T'][1] = -1;
  lookupTable[(unsigned char)'t'][0] = 3;
  lookupTable[(unsigned char)'t'][1] = -1;

  for(int i = 0; i < this->k; i++) {
    this->kmask = this->kmask << 2;
    this->kmask |= 3;
  }
  this->size = 0;
}

void Hash::intialize(string s) {
  this->pos = 0;
  this->len = s.length();
  this->seq = s;
  this->size = 0;
}

size_t Hash::nextHash(unsigned long long int &hashcode) {
  if(pos == len)
    return -1;

  while(1) {
    if(pos == len)
      return -1;
    hashcode = hashcode << 2;
    hashcode |= lookupTable[(unsigned char)this->seq[pos]][0];
    hashcode &= lookupTable[(unsigned char)this->seq[pos]][1];
    this->size++;
    this->size &= lookupTable[(unsigned char)this->seq[pos]][1];

    if(size >= k) {
      hashcode = hashcode & this->kmask;
      pos++;
      return 0;
    }
    pos++;
  }
}
