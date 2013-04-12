// enveomics/regex.h - Regular Expressions for enve-omics software
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @license artistic 2.0
// @version 1.0

#ifndef ENVEOMICS_REGEX_H
#define ENVEOMICS_REGEX_H

#include <regex.h>

using namespace std;

// The rreplace() function from: http://www.daniweb.com/software-development/c/code/216955
int rreplace (char *buf, int size, regex_t *re, char *rp);

// The str_replace() function from: http://www.zedwood.com/article/105/cpp-strreplace-function
string& str_replace(const string &search, const string &replace, string &subject);

#endif

