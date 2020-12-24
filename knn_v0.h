#ifndef knn_v0_h__
#define knn_v0_h__

#include "auxlib.h"

/**********************************************
 *            Debugging Options
 **********************************************/
#define DEBUG_PRINT_D1  false
#define DEBUG_PRINT_D2  false
#define DEBUG_PRINT_D3  false
/**********************************************/


/**
 * Reads m Query Points, with d Dimennsions in 
 * Row-Major Format, and finds the k nearest points
 * for each of those m points, in a set of n Corpus
 * Points X.
 **/
knnresult kNN(double * X, double * Y, int n, int m, int d, int k);


#endif

