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


#endif

