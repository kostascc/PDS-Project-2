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
#include <cilk/cilk.h>



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
 * True if V1 is running.
 * Default: False.
 **/
extern bool _MODE_V1_RUNNING ;

/**
 * Print kNN Neighbors result.
 * Default: False.
 **/
extern bool _KNN_PRINT ;

/**
 * Print calculated distance matrix.
 * Default: False.
 **/
extern bool _DIST_PRINT ;

/**
 * Print timing information.
 * Default: true.
 **/
extern bool _TIMER_PRINT ;



/**
 * kNN result struct
 **/
typedef struct knnresult{
  int    * nidx;    //!< Indices (0-based) of nearest neighbors [m-by-k]
  double * ndist;   //!< Distance of nearest neighbors          [m-by-k]
  int      m;       //!< Number of query points                 [scalar]
  int      k;       //!< Number of nearest neighbors            [scalar]
} knnresult;



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



#endif