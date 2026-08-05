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
#include <vector>
#include <cstring>

void init_multinode(int& argc, char**& argv, int& pid, int& pp) {
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &processID);
}
void finalize_multinode() {
  MPI_Finalize();
}

void barrier_multinode() {
  MPI_Barrier(MPI_COMM_WORLD);
}

void broadcast_bool(void* value) {
  MPI_Bcast(value, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  barrier_multinode();
}

void broadcast_int(void* value) {
  MPI_Bcast(value, 1, MPI_INT, 0, MPI_COMM_WORLD);
  barrier_multinode();
}

void broadcast_size_t(void* value) {
  MPI_Bcast(value, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  barrier_multinode();
}

void broadcast_double(void* value) {
  MPI_Bcast(value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  barrier_multinode();
}

void broadcast_char(void* value, size_t size) {
  MPI_Bcast(value, size, MPI_CHAR, 0, MPI_COMM_WORLD);
  barrier_multinode();
}
void broadcast_char(void* value) {
  MPI_Bcast(value, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
  barrier_multinode();
}

void reduce_sum_int(int *send, int *receive, int size){
  MPI_Reduce(send, receive, size, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
}
// BUG (not fixed here, just flagged): `receive` is passed by value, so the
// result computed below is discarded when the function returns -- the
// caller never sees it. Currently unused anywhere in this codebase (every
// call site uses the array overload above), so it's a latent bug rather
// than an active one, but worth fixing (change to `int &receive`) or
// removing before anything starts relying on it.
void reduce_sum_int(int send, int receive){
  int *send_ar = new int[1], *receive_ar = new int[1];
  send_ar[0] = send;
  MPI_Reduce(send_ar, receive_ar, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  receive = receive_ar[0];
}

void reduce_sum_double(double *send, double *receive, int size){
  MPI_Reduce(send, receive, size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

int ranks_on_this_node() {
  char name[MPI_MAX_PROCESSOR_NAME];
  int  len;
  MPI_Get_processor_name(name, &len);

  std::vector<char> all_names(
    (size_t) processes * MPI_MAX_PROCESSOR_NAME
  );
  MPI_Allgather(
    name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR,
    all_names.data(), MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MPI_COMM_WORLD
  );

  int count = 0;
  for (int i = 0; i < processes; i++)
    if (strncmp(
          &all_names[(size_t) i * MPI_MAX_PROCESSOR_NAME], name,
          MPI_MAX_PROCESSOR_NAME
        ) == 0) count++;

  return count;
}

#else
void init_multinode(int& argc, char**& argv, int& pid, int& pp){
  pid = 0;
  pp = 1;
}
void finalize_multinode() {}
void broadcast_bool(void* value) {}
void broadcast_int(void* value) {}
void broadcast_size_t(void* value) {}
void broadcast_double(void* value) {}
void broadcast_char(void* value, size_t size) {}
void broadcast_char(void* value) {}
void barrier_multinode() {}
void reduce_sum_int(int *send, int *receive, int size) {}
void reduce_sum_int(int send, int receive) {}
void reduce_sum_double(double *send, double *receive, int size) {}
int ranks_on_this_node() { return 1; }
#endif

