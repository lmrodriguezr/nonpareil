// go - Gene Ontology Library for enve-omics software
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @license artistic 2.0
// @version 1.0

#ifndef ENVEOMICS_GO_H
#define ENVEOMICS_GO_H

enum go2sqlite_option {
  GO2SQLITE_NOEXTDB = 0x01,
  GO2SQLITE_NOSYNONYM = 0x02,
  GO2SQLITE_NODEF = 0x04
};

enum godomain {
   GODOMAIN_BP = 1,
   GODOMAIN_MF = 2,
   GODOMAIN_CC = 3
};

/**
 * int go2sqlite(char *sqlFile, char *goDb, int opts);
 * Description:
 *   Takes a MySQL-based SQL file containing the GO database and creates an SQLite3 database.
 * Input:
 *   char *sqlFile: Path to the MySQL-based SQL with the GO DB.
 *   char *goDb: Path to the SQLite3 database to be created.
 *   int opts: Zero or more options.  Options are defined in the enum go2sqlite_option, and
 *      can be combined with | (the bit-wise "or").
 * Output:
 *   Returns zero on success, a non-zero integer otherwise.
 */
int go2sqlite(char *sqlFile, char *goDb, int opts);

/**
 * godomain go_get_domain(char *goDb, unsigned int go_id);
 * Description:
 *   Gets the GO Domain of a given GO ID.
 * Input:
 *   char *goDb: Path to the SQLite3 database, created by go2sqlite().
 *   size_t go_id: The GO ID.
 */
godomain go_get_domain(char *goDb, size_t go_id);

int go_length(char *goDb, size_t go_id_a, size_t go_id_b);
int go_get_nodes_parents(char *goDb, size_t *parents, size_t *go_ids, size_t nodes_no);
int go_get_node_parents(char *goDb, size_t *parents, size_t go_id);

#endif

