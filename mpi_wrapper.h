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
 * @brief Initializes MPI Communications
 * 
 * @param node_id 
 * @param cluster_size 
 */
void mpi_initialize(int* node_id, int* cluster_size);

/**
 * @brief Aborts MPI Execution in COMM_WORLD.
 * (The Whole job stops)
 **/
void mpi_abort();


/**
 * @brief Aborts MPI Execution in COMM_WORLD,
 * with an error message.
 * (The Whole job stops)
 * 
 * @param msg 
 */
void mpi_abort_msg(char* msg);

/**
 * @brief Finishes local job silently
 */
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
 * @brief Non-Blocking Send with DataType.
 * 
 * @param code 
 * @param Y 
 * @param m 
 * @param partner 
 * @param type 
 * @param request 
 */
void mpi_send_data_nb_t(int code, double* Y, int m, int partner, MPI_Datatype type, MPI_Request request[]);


/**
 * @brief Non-Blocking Receive with DataType.
 * 
 * @param code 
 * @param Z 
 * @param m 
 * @param partner 
 * @param type 
 * @param request 
 */
void mpi_receive_data_nb_t(int code, double* Z, int m, int partner, MPI_Datatype type, MPI_Request request[]);


/**
 * @brief Non-Blocking Send,
 * with default DataType: Double.
 * 
 * @param code 
 * @param Y 
 * @param m 
 * @param partner 
 * @param request 
 */
void mpi_send_data_nb(int code, double* Y, int m, int partner, MPI_Request request[]);


/**
 * @brief Non-Blocking Receive,
 * with default DataType: Double.
 * 
 * @param code 
 * @param Z 
 * @param m 
 * @param partner 
 * @param request 
 */
void mpi_receive_data_nb(int code, double* Z, int m, int partner, MPI_Request request[]);


/**
 * @brief Blocks until send has finished.
 * 
 * @param request 
 */
void mpi_send_data_wait(MPI_Request request[]);


/**
 * @brief Blocks until receive has finished.
 * 
 * @param request 
 */
void mpi_receive_data_wait(MPI_Request request[]);


/**
 * @brief Blocking Send,
 * with default DataType: Double.
 * 
 * @param code 
 * @param Y 
 * @param m 
 * @param partner 
 * @param request 
 */
void mpi_send_data_b(int code, double* Y, int m, int partner, MPI_Request request[]);


/**
 * @brief Blocking Receive,
 * with default DataType: Double.
 * 
 * @param code 
 * @param Z 
 * @param m 
 * @param partner 
 * @param request 
 */
void mpi_receive_data_b(int code, double* Z, int m, int partner, MPI_Request request[]);


/**
 * @brief Blocking Send with DataType.
 * 
 * @param code 
 * @param Y 
 * @param m 
 * @param partner 
 * @param type 
 * @param request 
 */
void mpi_send_data_b_t(int code, double* Y, int m, int partner, MPI_Datatype type, MPI_Request request[]);


/**
 * @brief Blocking Receive with DataType.
 * 
 * @param code 
 * @param Z 
 * @param m 
 * @param partner 
 * @param type 
 * @param request 
 */
void mpi_receive_data_b_t(int code, double* Z, int m, int partner, MPI_Datatype type, MPI_Request request[]);

/*********************************************/

#endif