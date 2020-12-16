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



knnresult knn_v1(double * X, double * Y, int n, int m, int d, int k);


#endif

