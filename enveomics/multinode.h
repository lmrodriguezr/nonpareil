// enveomics/mpi.h - Library for mpi-related code in enve-omics software
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @licanse artistic 2.0
// @version 1.0


#ifndef ENVEOMICS_MULTINODE_H
#define ENVEOMICS_MULTINODE_H

/**
 * void init_multinode(int argc, char **argv, int &processID, int &processes);
 * Description:
 *   Initializes multinode mode. If MPI is not supported, it's a dummy function.
 * Input:
 *   int argc: argc of main.
 *   char **argv: argv of main.
 *   int &pid: Reference to an integer to be defined as the (current) process ID.
 *   int &pp: Reference to an integer to be defined as the total number of spanned processes.
 */
void init_multinode(int& argc, char**& argv, int& pid, int& pp);

/**
 * void finalize_multinode();
 * Description:
 *   If MPI support is active closes MPI, otherwise it's a dummy function.
 */
void finalize_multinode();

size_t	broadcast_int(size_t value);
double	broadcast_double(double value);
char*	broadcast_char(char* value, size_t size);
char	broadcast_char(char value);
void	barrier_multinode();
void	reduce_sum_int(int *send, int *receive, int size);

#endif

