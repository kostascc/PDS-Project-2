#ifndef mpi_wrapper_h_
#define mpi_wrapper_h_


#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>


/**********************************************
 *            Debugging Options
 **********************************************/
#define DEBUG_MPI                       false
/**********************************************/


/**
 * Initialize MPI Communications
 **/
void mpi_initialize(int* node_id, int* cluster_size);

/**
 * Aborts MPI Execution in COMM_WORLD.
 * (The Whole job stops)
 **/
void mpi_abort();

/**
 * Aborts MPI Execution in COMM_WORLD,
 * with an error message.
 * (The Whole job stops)
 **/
void mpi_abort_msg(char* msg);

/**
 * Finishes local job silently.
 **/
void mpi_finish_local();


/*********************************************
 *               MPI Calls
 *********************************************/

/**
 * Usage:
 * 
 * For non blocking:
 *   mpi_send_data_nb(MPI_MODE_DISTRIBUTING, X, n, node_send, mpi_request);
 *   mpi_send_data_wait(mpi_request); // call later
 * 
 * For blocking:
 *   mpi_send_data_b(...)
 * 
 * Use wait at any time to block while receiving or sending:
 *   mpi_send_wait(mpi_request); // Wait for sending to finish
 * 
 **/

/**
 * Non-Blocking Send with DataType.
 **/
void mpi_send_data_nb_t(int code, double* Y, int m, int partner, MPI_Datatype type, MPI_Request request[]);

/**
 * Non-Blocking Receive with DataType.
 **/
void mpi_receive_data_nb_t(int code, double* Z, int m, int partner, MPI_Datatype type, MPI_Request request[]);

/**
 * Non-Blocking Send,
 * with default DataType: Double.
 **/
void mpi_send_data_nb(int code, double* Y, int m, int partner, MPI_Request request[]);

/**
 * Non-Blocking Receive,
 * with default DataType: Double.
 **/
void mpi_receive_data_nb(int code, double* Z, int m, int partner, MPI_Request request[]);

/**
 * Blocks until send has finished.
 **/
void mpi_send_data_wait(MPI_Request request[]);

/**
 * Blocks until receive has finished.
 **/
void mpi_receive_data_wait(MPI_Request request[]);

/**
 * Blocking Send,
 * with default DataType: Double.
 **/
void mpi_send_data_b(int code, double* Y, int m, int partner, MPI_Request request[]);

/**
 * Blocking Receive,
 * with default DataType: Double.
 **/
void mpi_receive_data_b(int code, double* Z, int m, int partner, MPI_Request request[]);

/**
 * Blocking Send with DataType.
 **/
void mpi_send_data_b_t(int code, double* Y, int m, int partner, MPI_Datatype type, MPI_Request request[]);

/**
 * Blocking Receive with DataType.
 **/
void mpi_receive_data_b_t(int code, double* Z, int m, int partner, MPI_Datatype type, MPI_Request request[]);

/*********************************************/

#endif