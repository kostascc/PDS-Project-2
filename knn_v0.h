#ifndef knn_v0_h__
#define knn_v0_h__

#include "auxlib.h"
#include "mat.h"
#include <cblas.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <cilk/cilk.h>



knnresult knn_v0(double * X, double * Y, int n, int m, int d, int k);


#endif

