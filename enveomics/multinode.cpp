// enveomics/multinode.h - Library for mpi-related code in enve-omics software
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @licanse artistic 2.0
// @version 1.0

#include "universal.h"
#include "multinode.h"

extern int processID;
extern int processes;

#ifdef ENVEOMICS_MULTI_NODE
#include <mpi.h>

void init_multinode(int& argc, char**& argv, int& pid, int& pp){
   MPI::Init(argc, argv);
   processes = MPI::COMM_WORLD.Get_size();
   processID = MPI::COMM_WORLD.Get_rank();
}
void finalize_multinode(){
   MPI::Finalize();
}

void barrier_multinode(){
   MPI::COMM_WORLD.Barrier();
}

size_t broadcast_int(size_t value){
   int buffer[1];
   buffer[0]=value;
   MPI::COMM_WORLD.Bcast(buffer, 1, MPI::INT, 0);
   barrier_multinode();
   value=buffer[0];
   return(value);
}

double broadcast_double(double value){
   double buffer[1];
   buffer[0]=value;
   MPI::COMM_WORLD.Bcast(buffer, 1, MPI::DOUBLE, 0);
   barrier_multinode();
   value=buffer[0];
   return(value);
}

char* broadcast_char(char* value, size_t size){
   MPI::COMM_WORLD.Bcast(value, size, MPI::CHAR, 0);
   barrier_multinode();
   return(value);
}
char broadcast_char(char value){
   char buffer[1] = {value};
   MPI::COMM_WORLD.Bcast(buffer, 1, MPI::CHAR, 0);
   barrier_multinode();
   value=buffer[0];
   return(value);
}

void reduce_sum_int(int *send, int *receive, int size){
   MPI::COMM_WORLD.Reduce(send, receive, size, MPI::INT, MPI::SUM, 0);
}
void reduce_sum_int(int send, int receive){
   int *send_ar = new int[1], *receive_ar = new int[1];
   send_ar[0] = send;
   MPI::COMM_WORLD.Reduce(send_ar, receive_ar, 1, MPI::INT, MPI::SUM, 0);
   receive = receive_ar[0];
}

void reduce_sum_double(double *send, double *receive, int size){
   MPI::COMM_WORLD.Reduce(send, receive, size, MPI::DOUBLE, MPI::SUM, 0);
}

#else
void init_multinode(int& argc, char**& argv, int& pid, int& pp){
   pid = 0;
   pp = 1;
}
void finalize_multinode(){}
size_t	broadcast_int(size_t value){return(value);}
double	broadcast_double(double value){return(value);}
char*	broadcast_char(char* value, size_t size){return(value);}
char	broadcast_char(char value){return(value);}
void	barrier_multinode(){}
void	reduce_sum_int(int *send, int *receive, int size){}
void	reduce_sum_int(int send, int receive){}
void	reduce_sum_double(double *send, double *receive, int size){}
#endif

