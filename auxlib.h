#ifndef auxlib_h__
#define auxlib_h__


#include <stdlib.h>
#include <cilk/cilk.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mat.h"
#include <cblas.h>
#include <cilk/cilk.h>
#include "mpi_wrapper.h"
#include <mpi.h>



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
#define MPI_MODE_QUERY_DISTRIBUTION     90
/**********************************************/



/**
 * Environment variable name for
 * Calculated distance matrix prints.
 **/
#define _DIST_PRINT_VAR "DIST_PRINT"

/**
 * Environment variable name for
 * kNN neighbor print.
 **/
#define _KNN_PRINT_VAR "KNN_PRINT"

/**
 * Environment variable name for
 * timing info.
 **/
#define _TIMER_PRINT_VAR "TIMER_PRINT"


#define min(x,y) (((x) < (y)) ? (x) : (y))

#define max(x,y) (((x) > (y)) ? (x) : (y))


/**
 * Startup Script Options
 **/
#define MATRIX_PRINT  false  // Print Resulting X & Y
#define RAND_SEED     false // Use seed for randomization
#define MAX_VECTOR    40.0  // Max distance vector


/**
 * Each external variable should
 * be defined both in .c and .h files.
 **/

/**
 * @brief True if V1 is running.
 * Default: False.
 **/
extern bool _MODE_V1_RUNNING ;

/**
 * @brief True if V2 is running.
 * Default: False.
 */
extern bool _MODE_V2_RUNNING ;

/**
 * @brief Print kNN Neighbors result.
 * Default: False.
 **/
extern bool _KNN_PRINT ;

/**
 * @brief Print calculated distance matrix.
 * Default: False.
 **/
extern bool _DIST_PRINT ;

/**
 * @brief Print timing information.
 * Default: true.
 **/
extern bool _TIMER_PRINT ;


/**
 * @brief kNN result struct
 * 
 * @param nidx Indices (0-based) of nearest neighbors [m-by-k]
 * @param ndist Distance of nearest neighbors [m-by-k]
 * @param m Number of query points
 * @param k Number of nearest neighbors
 */
typedef struct knnresult{
  int    * nidx;    //!< Indices (0-based) of nearest neighbors [m-by-k]
  double * ndist;   //!< Distance of nearest neighbors          [m-by-k]
  int      m;       //!< Number of query points                 [scalar]
  int      k;       //!< Number of nearest neighbors            [scalar]
} knnresult;



/**
 * @brief Main File runtime struct
 * 
 * @param X
 * @param Y
 * @param n
 * @param m
 * @param d
 * @param k
 **/
typedef struct _runtime{
  double* X;
  double* Y;
  int n;
  int m;
  int d;
  int k;
} _runtime;


/**
 * Sorts and selects neighbors
 * for a specific m.
 **/
int aux_sort_idx(double** C, int** nidx, double** ndist, int N, int M, int m, int k);


/** 
 * Makes a copy of an array from 
 * start to end - 1. Equivalent to Python's 
 * arr[start:end] 
 **/
// Double Array
double * aux_slice_d(double ** arr, int start, int end);
// Integer Array
int * aux_slice_i(int ** arr, int start, int end);


/** 
 * Merge left and right into result, 
 * overwriting the array in result. 
 **/
void aux_merge ( double** _I, int ** _J, double ** left_I, int ** left_J, double ** right_I, int ** right_J, int leftLen, int rightLen); 


/**
 * In-Place Merge-Sort, with
 * CilkPlus parallelization.
 **/
void aux_mergeSort( double ** _I, int ** _J, int len );

/**
 * Print knn result
 **/
void print_res(knnresult* knn);


/**
 * Startup Script
 **/
_runtime startup(int argc, char** argv);

// TODO: For Debugging, Remove later
extern int tmp_node_id;


/**
 * @brief For Worker nodes, it receives results from previous nodes,
 * compares them to the local KNNs, and propagates the
 * best results clockwise.
 * For the Master node, it receives results from the last
 * worker node and returnes the best results in 'fin'. 
 * Warning: The 'm' & 'k' in each knnresult should be
 * correct.
 * 
 * @param res KNN Result Array
 * @param[out] fin Final KNN result (m and k are required to be set!)
 * @param node_id Local Node ID
 * @param receive_node Node to receive from
 * @param result_count Results to receive (in this case the batch count)
 */
void collect_n_propagate_knn(knnresult* res, knnresult* fin, int node_id, int node_receive, int node_send, int result_count);


void knnresult_check_offset(knnresult* res, const int min, const int max);


void knnresult_check_batch(double *Y, double *X, int length, int offset);


void knnresult_offset_nidx(knnresult* res, int offset);


int n_per_node(int node_id, int cluster_size, int n);


void compare_knnresult(knnresult* res, knnresult* res_);


#endif