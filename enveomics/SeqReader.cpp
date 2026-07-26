#include "SeqReader.h"
#include "universal.h"

using namespace std;
//remove the global variable
const double PhredQual [43] = {
  1.0,
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
  0.000063095734448
};

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

void buildFastqSeq(
      string header, string sequence, string qual, Sequence &out,
      int phredOffset) {
  vector<double> baseProb;

  for(size_t i = 0; i < qual.length(); i++) {
    int idx = (int)(unsigned char) qual[i] - phredOffset;
    // Clamp defensively: an unexpected quality character (wrong encoding
    // guess, corrupted data) must not read past PhredQual's bounds.
    if (idx < 0) idx = 0;
    if (idx >= 43) idx = 42;
    baseProb.push_back(PhredQual[idx]);
  }
  Sequence seq = Sequence(header, sequence, qual, baseProb);
  out = seq;
}

int detect_phred_offset(ifstream &ifs) {
  std::streampos start = ifs.tellg();
  const size_t MAX_RECORDS = 10000;
  size_t records_seen = 0, nline = 0;
  int minChar = 127;
  string line;

  ifs.clear();
  ifs.seekg(0, std::ios::beg);
  while (records_seen < MAX_RECORDS && getline(ifs, line)) {
    if (nline % 4 == 3) {
      for (size_t i = 0; i < line.length(); i++) {
        int c = (int)(unsigned char) line[i];
        if (c < minChar) minChar = c;
      }
      records_seen++;
    }
    nline++;
  }
  ifs.clear();
  ifs.seekg(start);

  // Sanger/Illumina 1.8+ (Phred+33) reaches as low as '!' (33); Illumina
  // 1.3-1.7 (Phred+64) practically never goes below '@'-ish (64). A
  // low quality character below 64 can only occur in the +33 encoding.
  return (minChar < 64) ? 33 : 64;
}

void buildFastaSeq(string header, string sequence, Sequence &out) {
  Sequence seq = Sequence(header, sequence);
  out = seq;
}

SeqReader::SeqReader(ifstream &ifs_, unsigned int rseed) : ifs(ifs_) {
  if (rseed == 0) rseed = rd();
  this->readNext = true;
  this->gen = std::mt19937_64(rseed);
  unsigned long long int end;
  this->ifs.seekg(0,std::ios::end);
  end = this->ifs.tellg();
  this->distribution =
    std::uniform_int_distribution<unsigned long long int>(0, end);
  ifs.seekg(0,std::ios::beg);
}

void SeqReader::reset() {
  this->ifs.clear();
  this->ifs.seekg(0,std::ios::beg);
  this->readNext = true;
}

FastqReader::FastqReader(ifstream &ifs) : SeqReader(ifs) {
  this->phredOffset = detect_phred_offset(this->ifs);
}
FastaReader::FastaReader(ifstream &ifs) : SeqReader(ifs) {}
FastqReader::FastqReader(ifstream &ifs, unsigned int rseed) : SeqReader(ifs, rseed) {
  this->phredOffset = detect_phred_offset(this->ifs);
}
FastaReader::FastaReader(ifstream &ifs, unsigned int rseed) : SeqReader(ifs, rseed) {}

size_t FastqReader::readNextSeq(Sequence &out) {
  if (this->readNext == false) reset();

  string temp;
  string header;
  string sequence;
  string qual;

  if (getline(this->ifs, header).eof())
    return -1;
    // Be careful: this actually returns the largest unsigned integer,
    // not -1, since the function's return type is size_t

  char c;

  if (!getline(this->ifs, sequence).good())
    error("The file does not have proper fastq format: missing sequence");
  this->ifs.get(c);
  if (c != '+')
    error("The file does not have proper fastq format: wrong separator");
  if (!getline(this->ifs, temp).good())
    error("The file does not have proper fastq format: missing separator");
  if (!getline(this->ifs, qual).good())
    error("The file does not have proper fastq format: missing quality");

  buildFastqSeq(header, sequence, qual, out, this->phredOffset);

  return 0;
}

size_t FastaReader::readNextSeq(Sequence &out) {
  if (this->readNext == false) reset();

  string header;
  string sequence;
  string temp;
  char c;
  if (getline(this->ifs, header).eof())
    return -1;
    // Be careful: this actually returns the largest unsigned integer,
    // not -1, since the function's return type is size_t
  while (true) {
    if (getline(this->ifs, temp).good()) {
      sequence = sequence + temp;
      if (!this->ifs.get(c).good()){
        ifs.seekg(-1,std::ios::cur);
        break;
      }
      if (c == '>'){
        ifs.seekg(-1,std::ios::cur);
        break;
      }
      ifs.seekg(-1,std::ios::cur);
    } else {
      return -1;
    }
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
  while (true) {
    while (true) {
      if (!getline(this->ifs, temp).good()) {
        this->ifs.clear();
        this->ifs.seekg(0, std::ios::beg);
        break;
      }
      this->ifs.get(c);
      if (c == '@') break;
      this->ifs.seekg(-1, std::ios::cur);
    }

    getline(this->ifs, header);
    if (this->randomHeaders.count(header) == 1) continue;

    if (!getline(this->ifs, sequence).good()) continue;
    this->ifs.get(c);
    if (c != '+') continue;
    if (!getline(this->ifs, temp).good()) continue;
    if (!getline(this->ifs, qual).good()) continue;
    break;
  }

  buildFastqSeq(header, sequence, qual, out, this->phredOffset);

  this->randomHeaders.insert(header);
  return 0;
}

size_t FastaReader::getRandomSeq(Sequence &out) {
  string header;
  string sequence;
  string temp;
  unsigned long long int pos = this->distribution(this->gen);

  this->readNext = false;
  this->ifs.seekg(pos, std::ios::beg);

  bool inSequence = false;
  int  rounds = 0;

  while (true) {
    if (!getline(this->ifs, temp).good()) {
      // Roll back to the beginning of the file
      inSequence = false;
      sequence = (string)"";
      this->ifs.clear();
      this->ifs.seekg(0, std::ios::beg);
      if (rounds++ > 2) error("Cannot find FastA sequences in file");
      continue;
    }

    if (temp[0] == '>') {
      if (inSequence) {
        break;
      } else {
        header = temp.substr(1);
        if (this->randomHeaders.count(header) != 1) inSequence = true;
      }
    } else if(inSequence) {
      sequence = sequence + temp;
    }
  }

  buildFastaSeq(header, sequence, out);

  this->randomHeaders.insert(header);
  return 0;
}
