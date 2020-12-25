#ifndef knn_v2_h_
#define knn_v2_h_

#include "auxlib.h"
#include "mpi_wrapper.h"


/**
 * VPTree node struct
 **/
typedef struct _node{
  /**
   * TODO: 
   **/
} _node;

knnresult distrAllkNNVPT(double * X, int n_all, int d, int k);


#endif