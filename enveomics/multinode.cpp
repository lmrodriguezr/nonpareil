// enveomics/multinode.h - Library for mpi-related code in enve-omics software
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @licanse artistic 2.0
// @version 2.0

#include "universal.h"
#include "multinode.h"

extern int processID;
extern int processes;

#ifdef ENVEOMICS_MULTI_NODE
#include <mpi.h>

void init_multinode(int& argc, char**& argv, int& pid, int& pp){
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &processID);
}
void finalize_multinode(){
  MPI_Finalize();
}

void barrier_multinode(){
  MPI_Barrier(MPI_COMM_WORLD);
}

void broadcast_int(void* value){
  MPI_Bcast(value, 1, MPI_INT, 0, MPI_COMM_WORLD);
  barrier_multinode();
}

void broadcast_double(void* value){
  MPI_Bcast(value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  barrier_multinode();
}

void broadcast_char(void* value, size_t size){
  MPI_Bcast(value, size, MPI_CHAR, 0, MPI_COMM_WORLD);
  barrier_multinode();
}
void broadcast_char(void* value){
  MPI_Bcast(value, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
  barrier_multinode();
}

void reduce_sum_int(int *send, int *receive, int size){
  MPI_Reduce(send, receive, size, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
}
void reduce_sum_int(int send, int receive){
  int *send_ar = new int[1], *receive_ar = new int[1];
  send_ar[0] = send;
  MPI_Reduce(send_ar, receive_ar, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  receive = receive_ar[0];
}

void reduce_sum_double(double *send, double *receive, int size){
  MPI_Reduce(send, receive, size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

#else
void init_multinode(int& argc, char**& argv, int& pid, int& pp){
  pid = 0;
  pp = 1;
}
void finalize_multinode() {}
void broadcast_int(void* value) {}
void broadcast_double(void* value) {}
void broadcast_char(void* value, size_t size) {}
void broadcast_char(void* value) {}
void barrier_multinode() {}
void reduce_sum_int(int *send, int *receive, int size) {}
void reduce_sum_int(int send, int receive) {}
void reduce_sum_double(double *send, double *receive, int size) {}
#endif

