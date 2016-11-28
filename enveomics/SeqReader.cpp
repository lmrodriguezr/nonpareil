#include "SeqReader.h"
#include "universal.h"

using namespace std;
//remove the global variable
const double PhredQual [43] =
 {1.0,
  0.794328234724,
  0.63095734448,
  0.501187233627,
  0.398107170553,
  0.316227766017,
  0.251188643151,
  0.199526231497,
  0.158489319246,
  0.125892541179,
  0.1,
  0.0794328234724,
  0.063095734448,
  0.0501187233627,
  0.0398107170553,
  0.0316227766017,
  0.0251188643151,
  0.0199526231497,
  0.0158489319246,
  0.0125892541179,
  0.01,
  0.00794328234724,
  0.0063095734448,
  0.00501187233627,
  0.00398107170553,
  0.00316227766017,
  0.00251188643151,
  0.00199526231497,
  0.00158489319246,
  0.00125892541179,
  0.001,
  0.000794328234724,
  0.00063095734448,
  0.000501187233627,
  0.000398107170553,
  0.000316227766017,
  0.000251188643151,
  0.000199526231497,
  0.000158489319246,
  0.000125892541179,
  0.0001,
  0.0000794328234724,
  0.000063095734448};

std::random_device rd;

Sequence::Sequence(string header, string sequence, string qualstr, vector<double> baseProb) {
  this->header = header;
  this->sequence = sequence;
  this->qualstr = qualstr;
  this->baseProb = baseProb;
}

Sequence::Sequence(string header, string sequence) {
  this->header = header;
  this->sequence = sequence;
}

Sequence::Sequence() {
  this->header = "";
  this->sequence = "";
  this->qualstr = "";
}

void buildFastqSeq(string header, string sequence, string qual, Sequence &out) {
  vector<double> baseProb;

  for(size_t i = 0; i < qual.length(); i++) {
    baseProb.push_back(PhredQual[int(qual[i]) - 33]);
  }
  Sequence seq = Sequence(header, sequence, qual, baseProb);
  out = seq;
}

void buildFastaSeq(string header, string sequence, Sequence &out) {
  Sequence seq = Sequence(header, sequence);
  out = seq;
}

SeqReader::SeqReader(ifstream &ifs_) : ifs(ifs_) {
  this->readNext = true;
  this->gen = std::mt19937_64(rd());
  unsigned long long int end;
  this->ifs.seekg(0,std::ios::end);
  end = this->ifs.tellg();
  this->distribution = std::uniform_int_distribution<unsigned long long int>(0,end);
  ifs.seekg(0,std::ios::beg);
}

void SeqReader::reset() {
  this->ifs.clear();
  this->ifs.seekg(0,std::ios::beg);
  this->readNext = true;
}

FastqReader::FastqReader(ifstream &ifs) : SeqReader(ifs){}
FastaReader::FastaReader(ifstream &ifs) : SeqReader(ifs){}

size_t FastqReader::readNextSeq(Sequence &out) {
  if(this->readNext == false) {
    reset();
  }
  string temp;
  string header;
  string sequence;
  string qual;

  if(getline(this->ifs,header).eof())
    return -1;
    // Be careful: this actually returns the largest unsigned integer,
    // not -1, since the function's return type is size_t

  char c;
  if(getline(this->ifs,sequence).good()) {
    this->ifs.get(c);
    if(c != '+') error("The file provided does not have the proper fastq format");
  }
  else{
    error("The file provided does not have proper fastq format");
  }
  if(!getline(this->ifs,temp).good()) {
    error("The file you provided does not have the proper fastq format");
  }
  if(!getline(this->ifs,qual).good()) {
    error("The file you provided does not have the proper fastq format");
  }

  buildFastqSeq(header, sequence, qual,out);

  return 0;
}

size_t FastaReader:: readNextSeq(Sequence &out) {
  if(this->readNext == false) {
    reset();
  }

  string header;
  string sequence;

  if(getline(this->ifs,header).eof())
    return -1;
    // Be careful: this actually returns the largest unsigned integer,
    // not -1, since the function's return type is size_t

  if(!getline(this->ifs,sequence).good()) {
    error("This file provide does not have proper fasta format");
  }

  buildFastaSeq(header, sequence, out);

  return 0;
}

size_t FastqReader::getRandomSeq(Sequence &out) {
  string header;
  string sequence;
  string qual;

  unsigned long long int pos;
  pos = this->distribution(this->gen);
  this->readNext = false;
  this->ifs.seekg(pos,std::ios::beg);
  string temp;
  char c;
  while(true) {
    while(true) {
      if(!getline(this->ifs,temp).good()) {
        this->ifs.clear();
        this->ifs.seekg(0,std::ios::beg);
        break;
      }
      this->ifs.get(c);
      if(c == '@') {
        break;
      }
    }

    getline(this->ifs,header);
    if(this->randomHeaders.count(header) == 1)
      continue;
    getline(this->ifs,sequence);

    this->ifs.get(c);
    if(c != '+') {
      continue;
    }
    getline(this->ifs,temp);
    getline(this->ifs,qual);
    break;
  }

  buildFastqSeq(header,sequence,qual,out);

  this->randomHeaders.insert(header);
  return 0;
}

size_t FastaReader::getRandomSeq(Sequence &out) {
  string header;
  string sequence;

  unsigned long long int pos;
  pos = this->distribution(this->gen);
  this->readNext = false;
  this->ifs.seekg(pos,std::ios::beg);
  string temp;
  char c;
  while(true) {
    while(true) {
      if(!getline(this->ifs,temp).good()) {
        this->ifs.clear();
        this->ifs.seekg(0,std::ios::beg);
        break;
      }
      this->ifs.get(c);
      if(c == '>') {
        break;
      }
    }
    getline(this->ifs, header);
    if(this->randomHeaders.count(header) == 1)
      continue;
    getline(this->ifs,sequence);
    break;
  }

  buildFastaSeq(header,sequence,out);

  this->randomHeaders.insert(header);
  return 0;
}
