#ifndef auxlib_h__
#define auxlib_h__


#include <stdlib.h>
#include <cilk/cilk.h>



/**
 * kNN result struct
 **/
typedef struct knnresult{
  int    * nidx;    //!< Indices (0-based) of nearest neighbors [m-by-k]
  double * ndist;   //!< Distance of nearest neighbors          [m-by-k]
  int      m;       //!< Number of query points                 [scalar]
  int      k;       //!< Number of nearest neighbors            [scalar]
} knnresult;


#define min(x,y) (((x) < (y)) ? (x) : (y))

#define max(x,y) (((x) > (y)) ? (x) : (y))


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