// enveomics/sqlite.h - Library for SQLite usage in enve-omics software
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @license artistic 2.0
// @version 1.0

#ifndef ENVEOMICS_SQLITE_H
#define ENVEOMICS_SQLITE_H
#define SQLITE_MAX_COLS 400

#include <string>

using namespace std;

static int sqlite_callback(void *NotUsed, int argc, char **argv, char **azColName);
int sqlite_execute(char *dbFile, char *sql);
int sqlite_execute(char *dbFile, string sql);

#endif

