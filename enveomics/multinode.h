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

/**
 * void broadcast_bool(void* value);
 * void broadcast_int(void* value);
 * void broadcast_size_t(void* value);
 * void broadcast_double(void* value);
 *
 * Description:
 *   Broadcasts a single value of the corresponding type from rank 0 to
 *   every other rank, then synchronizes all ranks on a barrier. If MPI
 *   support is not active, these are dummy functions (no-ops). All ranks
 *   must call the same broadcast function, in the same order, for this to
 *   behave correctly.
 * Input:
 *   - `void *value`: Pointer to the value to broadcast. On rank 0, this is
 *     read; on every other rank, this is overwritten with rank 0's value
 */
void broadcast_bool(void* value);
void broadcast_int(void* value);
void broadcast_size_t(void* value);
void broadcast_double(void* value);

/**
 * void broadcast_char(void* value, size_t size);
 * void broadcast_char(void* value);
 *
 * Description:
 *   Broadcasts a buffer of `char` from rank 0 to every other rank (e.g. a
 *   fixed-size path buffer), then synchronizes all ranks on a barrier. The
 *   single-argument overload broadcasts exactly one `char`. If MPI support
 *   is not active, these are dummy functions (no-ops).
 * Input:
 *   - `void *value`: Pointer to the buffer to broadcast. On rank 0, this is
 *     read; on every other rank, this is overwritten with rank 0's buffer
 *   - `size_t size` (first overload only): Number of `char` elements in
 *     the buffer
 */
void broadcast_char(void* value, size_t size);
void broadcast_char(void* value);

/**
 * void barrier_multinode();
 *
 * Description:
 *   Blocks the calling rank until every rank has reached this point. If
 *   MPI support is not active, this is a dummy function (no-op).
 */
void barrier_multinode();

/**
 * void reduce_sum_int(int *send, int *receive, int size);
 * void reduce_sum_int(int send, int receive);
 *
 * Description:
 *   Sums a value (or array of values) element-wise across every rank, and
 *   stores the result on rank 0 only (other ranks' `receive` is left
 *   unset). If MPI support is not active, these are dummy functions
 *   (no-ops) -- do not rely on `receive` being populated when built
 *   without MPI. Every rank must call this with the same `size`.
 * Input:
 *   - `int *send` / `int send`: This rank's value(s) to contribute to the
 *     sum
 *   - `int *receive` / `int receive`: Where to store the summed result, on
 *     rank 0
 *   - `int size` (array overload only): Number of elements in `send`/
 *     `receive`
 */
void reduce_sum_int(int *send, int *receive, int size);
void reduce_sum_int(int send, int receive);

/**
 * void reduce_sum_double(double *send, double *receive, int size);
 *
 * Description:
 *   Sums an array of `double` element-wise across every rank, and stores
 *   the result on rank 0 only (other ranks' `receive` is left unset). If
 *   MPI support is not active, this is a dummy function (no-op). Every
 *   rank must call this with the same `size`.
 * Input:
 *   - `double *send`: This rank's values to contribute to the sum
 *   - `double *receive`: Where to store the summed result, on rank 0
 *   - `int size`: Number of elements in `send`/`receive`
 */
void reduce_sum_double(double *send, double *receive, int size);

/**
 * int ranks_on_this_node();
 * Description:
 *   Counts how many MPI ranks share the current rank's hostname. Outside MPI
 *   support, always returns 1.
 * Output:
 *   `int`: Number of ranks (including this one) running on the same node.
 */
int ranks_on_this_node();

#endif

