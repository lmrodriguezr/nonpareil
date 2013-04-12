// go.h - Gene Ontology Library for enve-omics software
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @license artistic 2.0
// @version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <regex.h>
#include <sqlite3.h>

#include "universal.h"
#include "regex.h"
#include "sqlite.h"
#include "go.h"

using namespace std;

#define DEBUG(a) (cerr << "(LINE " << a << ")" << endl)
// #define DEBUG(a) (a)
#define MAX_GO_NODES 2048
#define MAX_GO_SIBLINGS 128
#define MAX_GO_DISTANCE 128

int go2sqlite(char *sqlFile, char *goDb, int opts){
   // Vars
   regex_t rxCut, rxInsert, rxInsertDb, rxInsertSyn, rxInsertDef, rxValues, rx0, rx1, rx2, rx3, rx4, rx5;
   regmatch_t match;
   ifstream myfile;
   string line, fullSql, tmpstr;
   char *tmpchrarr;
   int executed=0;
   // Load SQL in memory
   myfile.open(sqlFile, ios::in);
   if(!myfile.is_open()){
      error("Impossible to open the file", sqlFile);
      return 1;
   }
   fullSql = "";
   regcomp(&rxCut, (char *)";",REG_EXTENDED);
   regcomp(&rxInsert, (char *)"INSERT INTO `[A-Za-z0-9_]*` VALUES ?",REG_EXTENDED);
   regcomp(&rxInsertDb, (char *)"INSERT INTO `(term_)?db(xref)?` VALUES ?",REG_EXTENDED);
   regcomp(&rxInsertSyn, (char *)"INSERT INTO `term_synonym` VALUES ?",REG_EXTENDED);
   regcomp(&rxInsertDef, (char *)"INSERT INTO `term_definition` VALUES ?",REG_EXTENDED);
   regcomp(&rxValues, (char *)"\\(([0-9\\.eE-]+|NULL|'([^']|'')+'|,)+\\)[,;]",REG_EXTENDED);
   regcomp(&rx0,(char *)"AUTO_INCREMENT=[0-9]+",REG_EXTENDED);
   regcomp(&rx1,(char *)"AUTO_INCREMENT",REG_EXTENDED);
   regcomp(&rx2,(char *)",\n * KEY [^\n]*)",REG_EXTENDED);
   regcomp(&rx3,(char *)",\n * UNIQUE KEY [^\n]*)",REG_EXTENDED);
   regcomp(&rx4,(char *)" *TYPE=MyISAM *;",REG_EXTENDED);
   //regcomp(&rx4,(char *)" *ENGINE=MyISAM *;",REG_EXTENDED);
   regcomp(&rx5,(char *)"U?N?LOCK TABLES[^;]*;",REG_EXTENDED);
   while(myfile.good()){
      getline(myfile, line);
      fullSql.append(line);
      fullSql.append((char *)"\n");
      // lrr Sep.11.11:
      // Unknown error caused the fullSql to go blank after copying into tmpchrarr.
      // The workaround: to copy fullSql into a new var (tmpchrarrI) and then copy
      // that new var into tmpchrarr.
      // {
      char tmpchrarrI [fullSql.size()];
      for(int a=0; a<=fullSql.size(); a++) tmpchrarrI[a] = fullSql[a];
      tmpchrarr = (char *) tmpchrarrI;
      // }
      if(regexec(&rxCut,tmpchrarr,1,&match,0)==0){
	 // Fix SQL to match SQLite requirements
	 if (  rreplace (tmpchrarr, fullSql.size(), &rx0, (char *)"") |
	       rreplace (tmpchrarr, fullSql.size(), &rx1, (char *)"") |
	       rreplace (tmpchrarr, fullSql.size(), &rx2, (char *)"") |
	       rreplace (tmpchrarr, fullSql.size(), &rx3, (char *)"") |
	       rreplace (tmpchrarr, fullSql.size(), &rx4, (char *)";") |
	       rreplace (tmpchrarr, fullSql.size(), &rx5, (char *)"") |
	       false ){
	    error("Regular expression failed");
	    return 1;
	 }
	 // And execute
	 if(regexec(&rxInsert, tmpchrarr,1,&match,0)==0){
	    int insertLen = match.rm_eo - match.rm_so;
	    char insert[insertLen], *values, *init_insert, *close_delete;
	    bool knownTable;
	    memmove(insert, tmpchrarr + match.rm_so, insertLen);
	    insert[insertLen] = NULL;
	    // Do I need to run this INSERT?
	    if(   (!( (GO2SQLITE_NOEXTDB && opts) & regexec(&rxInsertDb, tmpchrarr, 1, &match, 0)==0 & match.rm_so>=0)) &
	          (!( (GO2SQLITE_NOSYNONYM && opts) & regexec(&rxInsertSyn, tmpchrarr, 1, &match, 0)==0 & match.rm_so>=0)) &
	          (!( (GO2SQLITE_NODEF && opts) & regexec(&rxInsertDef, tmpchrarr, 1, &match, 0)==0 & match.rm_so>=0)) ){
	       // Type of INSERT
	       knownTable = true;
	       if(strcmp(insert, "INSERT INTO `graph_path` VALUES ")==0){
		  init_insert =
		     (char*)"INSERT INTO `graph_path` SELECT -1 AS `id`,1 AS `term1_id`,1 AS `term2_id`,1 AS `relationship_type_id`,1 AS `distance`,1 AS `relation_distance`\n";
		  close_delete = (char*)"; \nDELETE FROM `graph_path` WHERE `id`=-1;";
	       }else if(strcmp(insert, "INSERT INTO `term` VALUES ")==0){
		  init_insert = (char*)"INSERT INTO `term` SELECT -1 AS `id`,1 AS `name`,1 AS `term_type`,1 AS `acc`,1 AS `is_obsolete`,0 AS `is_root`,0 AS `is_relation`\n";
		  close_delete = (char*)"; \nDELETE FROM `term` WHERE `id`=-1;";
	       }else if(strcmp(insert, "INSERT INTO `term2term` VALUES ")==0){
		  init_insert = (char*)"INSERT INTO `term2term` SELECT -1 AS `id`,1 AS `relationship_type_id`,1 AS `term1_id`,1 AS `term2_id`,0 AS `complete`\n";
		  close_delete = (char*)"; \nDELETE FROM `term2term` WHERE `id`=-1;";
	       }else if(strcmp(insert, "INSERT INTO `term2term_metadata` VALUES ")==0){
		  init_insert = (char*)"INSERT INTO `term2term_metadata` SELECT -1 AS `id`,1 AS `relationship_type_id`,1 AS `term1_id`,1 AS `term2_id`\n";
		  close_delete = (char*)"; \nDELETE FROM `term2term_metadata` WHERE `id`=-1;";
	       }else if(strcmp(insert, "INSERT INTO `term_subset` VALUES ")==0){
		  init_insert = (char*)"INSERT INTO `term_subset` SELECT -1 AS `term_id`,-1 AS `subset_id`\n";
		  close_delete = (char*)"; \nDELETE FROM `term_subset` WHERE `term_id`=-1 AND `subset_id`=-1;";
	       }else{
		  knownTable = false;
	       }
	       
	       //if(strcmp(insert, "INSERT INTO `graph_path` VALUES ")==0){
	       if(knownTable){
	          // INSERT for known tables
		  int val_len = strlen(tmpchrarr), value_from=-1, insert_cols=0;
		  char in_string=0;
		  fullSql = (string)"";
		  fullSql.append(init_insert);
		  for(int pos = 0; pos<val_len; pos++){
		     if(in_string!=0){
		        if((pos==0 | tmpchrarr[pos-1]!='\\') & tmpchrarr[pos]==in_string) in_string = 0;
		     }else if((pos==0 | tmpchrarr[pos-1]!='\\') & (tmpchrarr[pos]=='\'' | tmpchrarr[pos]=='"')){
		        in_string = tmpchrarr[pos];
		     }else if(value_from>=0){
		        if(tmpchrarr[pos]==')'){
			   insert_cols++;
			   char one_row [pos-value_from];
			   memmove(one_row, tmpchrarr+value_from, pos-value_from);
			   one_row[pos-value_from] = NULL;
			   if(insert_cols >= SQLITE_MAX_COLS){
			      fullSql.append(close_delete);
			      fullSql.append(init_insert);
			      insert_cols = 0;
			      executed+=2;
			   }
			   fullSql.append(" UNION SELECT ");
			   fullSql.append(one_row);
			   value_from = -1;
			}
		     }else if(tmpchrarr[pos]=='('){
		        value_from = pos+1;
		     }

		     if(pos>0 & tmpchrarr[pos]=='\'' & tmpchrarr[pos-1]=='\\') tmpchrarr[pos-1]='\'';
		  }
		  fullSql.append(close_delete);
		  sqlite_execute(goDb, fullSql);
		  executed+=2;
		  cerr << " SQL statement " << executed << ": (modified) " << insert << "          \r";
	       }else{
	          // INSERT others
		  string tmpstrI = (string)tmpchrarr;
		  str_replace("\\'", "''", tmpstrI);
		  char valuesI [tmpstrI.size()];
		  for(int a=0; a<=tmpstrI.size(); a++) valuesI[a] = tmpstrI[a];
		  tmpchrarr = (char *)valuesI;
		  
		  while(regexec(&rxValues, tmpchrarr, 1, &match, 0)==0 && match.rm_so>=0){
		     values = new char [insertLen + match.rm_eo - match.rm_so];
		     memmove(values, insert, insertLen);
		     memmove(values+insertLen, tmpchrarr+match.rm_so, match.rm_eo - match.rm_so - 1);
		     values[insertLen + match.rm_eo - match.rm_so - 1] = NULL;
		     
		     for(int a=0; a<match.rm_eo; a++) tmpchrarr[a] = ' ';
		     sqlite_execute(goDb, values);
		     executed++;
		     if(executed % 10 == 0) cerr << " SQL statement " << executed << ": " << insert << "            \r";
		  }
	       }
	    }
	 }else{
	    // Otherwise just run it
	    sqlite_execute(goDb, tmpchrarr);
	    executed++;
	    if(executed % 10 == 0) cerr << " SQL statement " << executed << "                                               " << "\r";
	 }
	 fullSql = (string)"";
      }
   }
   myfile.close();
   regfree (&rxCut); regfree (&rxInsert); regfree (&rxValues);
   regfree (&rxInsertDb); regfree (&rxInsertSyn); regfree (&rxInsertDef);
   regfree (&rx0); regfree (&rx1); regfree (&rx2);
   regfree (&rx3); regfree (&rx4); regfree (&rx5);
   cerr << endl;
   return 0;
}

godomain go_get_domain(char *goDb, size_t go_id){
   sqlite3 *db;
   sqlite3_stmt *statement;
   char qry[100];
   char *str_domain;
   sprintf(qry, "SELECT `term_type` FROM `term` WHERE `acc`='GO:%07d';", (int)go_id);

   if(sqlite3_open(goDb, &db) != SQLITE_OK) error("Cannot open database", sqlite3_errmsg(db));
   if(sqlite3_prepare_v2(db, qry, -1, &statement, 0) != SQLITE_OK) error("Cannot prepare query", sqlite3_errmsg(db));
   if(sqlite3_step(statement) != SQLITE_ROW) error("Cannot find the GO ID", (int)go_id);
   str_domain = (char *)sqlite3_column_text(statement, 0);
   sqlite3_close(db);
   
   if(strcmp(str_domain, "molecular_function")==0) return GODOMAIN_MF;
   if(strcmp(str_domain, "biological_process")==0) return GODOMAIN_BP;
   if(strcmp(str_domain, "cellular_component")==0) return GODOMAIN_CC;
   error("Cannot find the GO Domain", str_domain);
}

int go_length(char *goDb, size_t go_id_a, size_t go_id_b){
   size_t	par_a[MAX_GO_SIBLINGS], par_b[MAX_GO_SIBLINGS],
   		par_a_t[MAX_GO_SIBLINGS], par_b_t[MAX_GO_SIBLINGS],
		nod_a[MAX_GO_SIBLINGS*MAX_GO_DISTANCE], nod_b[MAX_GO_SIBLINGS*MAX_GO_DISTANCE],
		nod_a_no=0, nod_b_no=0, par_a_no, par_b_no, length=0;
   
   if(go_id_a == go_id_b) return 0;
   
   par_a[0] = go_id_a;
   par_b[0] = go_id_b;

   while(true){
      if(par_a_no==0 || par_b_no==0) error("Impossible to find further ancestors, are these nodes in the same domain?", (double)go_id_b*1e-7+go_id_a);
      
      // Copy parents in the complete list of nodes (nod_*)
      // ToDo: I should avoid duplicates in nod_a and nod_b
      for(size_t i=0; i<par_a_no; i++) nod_a[nod_a_no++] = par_a[i];
      for(size_t i=0; i<par_b_no; i++) nod_b[nod_b_no++] = par_b[i];
      
      // Check if nod_a contains any of nod_b nodes
      for(size_t i=0; i<nod_a_no; i++)
         for(size_t j=0; i<nod_b_no; j++)
	    if(nod_a[i]==nod_b[j]){
	       free(par_a);
	       free(par_b);
	       free(par_a_t);
	       free(par_b_t);
	       free(nod_a);
	       free(nod_b);
	       return length;
	    }
      
      // Copy the parents (par_*) in temporal arrays (par_*_t)
      for(size_t i=0; i<par_a_no; i++) par_a_t[i] = par_a[i];
      for(size_t i=0; i<par_b_no; i++) par_b_t[i] = par_b[i];

      // Get the parents (in par_*) of the temporal arrays (par_*_t)
      par_a_no = go_get_nodes_parents(goDb, par_a, par_a_t, par_a_no);
      par_b_no = go_get_nodes_parents(goDb, par_b, par_b_t, par_b_no);
      
      // Increase the length
      length++;
      if(length>=MAX_GO_DISTANCE) error("Maximum distance between nodes reached", (int)length);
   }
}

int go_depth(char *goDb, size_t go_id){
   sqlite3 *db;
   sqlite3_stmt *statement;
   char qry[500];
   int	depth;

   sprintf(qry,
	"select g.distance from graph_path g, term t, term t2 where t.acc='GO:%07d' and t.id=g.term2_id and relationship_type_id=1 and t2.id=g.term1_id and t2.is_root=1;",
   	(int)go_id);

   if(sqlite3_open(goDb, &db) != SQLITE_OK) error("Cannot open database", sqlite3_errmsg(db));
   if(sqlite3_prepare_v2(db, qry, -1, &statement, 0) != SQLITE_OK) error("Cannot prepare query", sqlite3_errmsg(db));
   if(sqlite3_step(statement) != SQLITE_ROW) error("Cannot find the GO ID", (int)go_id);
   depth = sqlite3_column_int(statement, 0);
   sqlite3_close(db);

   return depth;
}

int go_get_nodes_parents(char *goDb, size_t *parents, size_t *go_ids, size_t nodes_no){
   size_t	*par_t, par_no=0, par_no_t;
   parents = (size_t *)malloc(((size_t)MAX_GO_NODES)*(sizeof(int)));
   
   for(size_t i=0; i<nodes_no; i++){
      par_no_t = go_get_node_parents(goDb, par_t, go_ids[i]);
      for(size_t j=0; j<par_no_t; j++) parents[par_no++] = par_t[i];
   }
   
   free(par_t);
   return par_no;
}

int go_get_node_parents(char *goDb, size_t *parents, size_t go_id){
   sqlite3 *db;
   sqlite3_stmt *statement;
   char 	qry[500];
   size_t	par_no=0;

   parents = (size_t *)malloc(((size_t)MAX_GO_SIBLINGS)*(sizeof(int)));
   sprintf(qry,
	"select t2.acc from graph_path g, term t, term t2 where t.acc='GO:%07d' and t.id=g.term2_id and relationship_type_id=1 and g.distance=1 and t2.id=g.term1_id;",
   	(int)go_id);

   if(sqlite3_open(goDb, &db) != SQLITE_OK) error("Cannot open database", sqlite3_errmsg(db));
   if(sqlite3_prepare_v2(db, qry, -1, &statement, 0) != SQLITE_OK) error("Cannot prepare query", sqlite3_errmsg(db));
   while(sqlite3_step(statement) == SQLITE_ROW)
      parents[par_no++] = sqlite3_column_int(statement, 0);
   sqlite3_close(db);

   return par_no;
}

