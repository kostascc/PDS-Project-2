#ifndef knn_v1_h__
#define knn_v1_h__

#include "auxlib.h"
#include "mat.h"
#include <cblas.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <cilk/cilk.h>
#include <mpi.h>
#include "knn_v0.h"


// Holds node ID,
// only for Debugging.
int _v1_node_id_local ;



/**********************************************
 *            Debugging Options
 **********************************************/
#define DEBUG_A             false
#define DEBUG_MPI           true
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


knnresult distrAllkNN(double * X, int n, int d, int k);


#endif

