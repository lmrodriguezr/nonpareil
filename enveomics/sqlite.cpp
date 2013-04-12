// enveomics/sqlite.h - Library for SQLite usage in enve-omics software
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @license artistic 2.0
// @version 1.0

#include <iostream>
#include <string>
#include <unistd.h>
#include <stdio.h>
#include <sqlite3.h>

#include "universal.h"

using namespace std;

static int sqlite_callback(void *NotUsed, int argc, char **argv, char **azColName){
   int i;
   for(i=0; i<argc; i++){
      cerr << azColName[i] << " = " << (argv[i] ? argv[i] : "NULL") << endl;
   }
   cerr << endl;
   return 0;
}

int sqlite_execute(char *dbFile, char *sql){
   sqlite3 *db;
   char *zErrMsg=0;
   if(sqlite3_open(dbFile, &db)!=SQLITE_OK) error("Can not open database", sqlite3_errmsg(db));
   if(sqlite3_exec(db, sql, sqlite_callback, 0, &zErrMsg)!=SQLITE_OK){
      sqlite3_close(db);
      error(zErrMsg, sql);
      return 1;
   }
   sqlite3_close(db);
   return 0;
}

int sqlite_execute(char *dbFile, string sql){
   char sqlChr[sql.size()];
   for(int a=0; a<=sql.size(); a++) sqlChr[a] = sql[a];
   return sqlite_execute(dbFile, sqlChr);
}

