#ifndef knn_v1_h__
#define knn_v1_h__

#include "auxlib.h"
#include "mpi_wrapper.h"
#include "knn_v0.h"


/**********************************************
 *            Debugging Options
 **********************************************/
#define DEBUG_A                         false
#define DEBUG_RES_I                     false
#define DEBUG_CHK_BATCH                 false
#define DEBUG_CHK_OFFSET                false
/**********************************************/


/**********************************************
 *          MPI Communication Tags
 **********************************************/
#define MPI_MODE_INITIALIZING           10
#define MPI_MODE_DATA_DISTRIBUTION      20
#define MPI_MODE_KNN_COLLECTION_DIST    30
#define MPI_MODE_KNN_COLLECTION_INDX    40
#define MPI_MODE_KNN_COLLECTION_M       50
#define MPI_MODE_CORPUS_DISTRIBUTION    60
#define MPI_MODE_CALCULATING            70
#define MPI_MODE_COLLECTING             80
/**********************************************/



/**
 * Reads n Query Points, with d Dimennsions in 
 * Row-Major Format, and finds the k nearest points
 * for each of those n points.
 **/
knnresult distrAllkNN(double * X, int n, int d, int k);


// /**
//  * Checks for indexes (can be used to check
//  * correct index offsetting based
//  * on node i).
//  **/
// void knnresult_check_offset(knnresult res, const int min, const int max);

// /**
//  * Checks Received Query Batch based on the
//  * whole corpus.
//  * length:'elements times dimensions'.
//  * Offset: the correct point in which the
//  * query can be found based on batch ID.
//  **/
// void knnresult_check_batch(double *Y, double *X, int length, int offset);

// /**
//  * Offsets indexes based on node ID.
//  **/
// void knnresult_offset_nidx(knnresult* res, int offset);

// /**
//  * Returns the number of corpus elements that
//  * shouuld be saved in a node (node ID), thus 
//  * equally the number of query points 
//  * per batch (batch ID).
//  **/
// int n_per_node(int node_id, int cluster_size, int n);

/**
 * Compare the results of different nodes on the same
 * Query, and find the resulting knn.
 * res  : Resulting data.
 * res_ : Comparison data.
 **/
// void compare_knnresult(knnresult* res, knnresult* res_);


#endif

