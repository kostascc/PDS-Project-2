#ifndef knn_v1_h__
#define knn_v1_h__

#include "auxlib.h"
#include <mpi.h>
#include "knn_v0.h"


/**********************************************
 *            Debugging Options
 **********************************************/
#define DEBUG_A             false
#define DEBUG_MPI           false
#define DEBUG_RES_I         false
#define DEBUG_CHK_BATCH     true
#define DEBUG_CHK_OFFSET    true
/**********************************************/


/**********************************************
 *          MPI Communication Tags
 **********************************************/
#define MPI_MODE_INITIALIZING           0
#define MPI_MODE_DATA_DISTRIBUTION      1
#define MPI_MODE_KNN_COLLECTION_DIST    2
#define MPI_MODE_KNN_COLLECTION_INDX    3
#define MPI_MODE_KNN_COLLECTION_M       7
#define MPI_MODE_CORPUS_DISTRIBUTION    6
#define MPI_MODE_CALCULATING            4
#define MPI_MODE_COLLECTING             5
/**********************************************/


/**
 * V1-Global Node ID,
 * Used for debugging only!
 **/
int _v1_node_id_local ;


/**
 * Reads n Query Points, with d Dimennsions in 
 * Row-Major Format, and finds the k nearest points
 * for each of those n points.
 **/
knnresult distrAllkNN(double * X, int n, int d, int k);


/**
 * Aborts MPI Execution in COMM_WORLD.
 * (The Whole job stops)
 **/
void abort();

/**
 * Finishes local job silently.
 **/
void finish_local();

/**
 * Checks for indexes (can be used to check
 * correct index offsetting based
 * on node i).
 **/
void _v1_check_offset(knnresult res, const int min, const int max);

/**
 * Checks Received Query Batch based on the
 * whole corpus.
 * m should be 'elements times dimensions'.
 * Offset is the correct point is which the
 * query can be found based on batch ID.
 **/
void _v1_check_batch(double *Y, double *X, int m, int offset);

/**
 * Offsets indexes based on node ID.
 **/
void _v1_offset_nidx(knnresult* res, int offset);

/**
 * Returns the number of corpus elements that
 * shouuld be saved in a node (node ID), thus 
 * equally the number of query points 
 * per batch (batch ID).
 **/
int _v1_n_per_node(int node_id, int cluster_size, int n);

/**
 * Compare the results of different nodes on the same
 * Query, and find the resulting knn.
 * res  : Resulting data.
 * res_ : Comparison data.
 **/
void _v1_compare_knnresult(knnresult* res, knnresult* res_);



/*********************************************
 *               MPI Calls
 *********************************************/

/**
 * Usage:
 * 
 * For non blocking:
 *   _v1_send_data_nb(MPI_MODE_DISTRIBUTING, X, n, node_send, mpi_request);
 *   _v1_send_data_wait(mpi_request); // call later
 * 
 * For blocking:
 *   _v1_send_data_b(...)
 * 
 * Use wait at any time to block while receiving or sending:
 *   _v1_send_wait(mpi_request); // Wait for sending to finish
 * 
 **/

/**
 * Non-Blocking Send with DataType.
 **/
void _v1_send_data_nb_t(int code, double* Y, int m, int partner, MPI_Datatype type, MPI_Request request[]);

/**
 * Non-Blocking Receive with DataType.
 **/
void _v1_receive_data_nb_t(int code, double* Z, int m, int partner, MPI_Datatype type, MPI_Request request[]);

/**
 * Non-Blocking Send,
 * with default DataType: Double.
 **/
void _v1_send_data_nb(int code, double* Y, int m, int partner, MPI_Request request[]);

/**
 * Non-Blocking Receive,
 * with default DataType: Double.
 **/
void _v1_receive_data_nb(int code, double* Z, int m, int partner, MPI_Request request[]);

/**
 * Blocks until send has finished.
 **/
void _v1_send_data_wait(MPI_Request request[]);

/**
 * Blocks until receive has finished.
 **/
void _v1_receive_data_wait(MPI_Request request[]);

/**
 * Blocking Send,
 * with default DataType: Double.
 **/
void _v1_send_data_b(int code, double* Y, int m, int partner, MPI_Request request[]);

/**
 * Blocking Receive,
 * with default DataType: Double.
 **/
void _v1_receive_data_b(int code, double* Z, int m, int partner, MPI_Request request[]);

/**
 * Blocking Send with DataType.
 **/
void _v1_send_data_b_t(int code, double* Y, int m, int partner, MPI_Datatype type, MPI_Request request[]);

/**
 * Blocking Receive with DataType.
 **/
void _v1_receive_data_b_t(int code, double* Z, int m, int partner, MPI_Datatype type, MPI_Request request[]);

/*********************************************/


#endif

