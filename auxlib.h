#ifndef auxlib_h__
#define auxlib_h__


// Definition of the kNN result struct
typedef struct knnresult{
  int    * nidx;    //!< Indices (0-based) of nearest neighbors [m-by-k]
  double * ndist;   //!< Distance of nearest neighbors          [m-by-k]
  int      m;       //!< Number of query points                 [scalar]
  int      k;       //!< Number of nearest neighbors            [scalar]
} knnresult;


#define min(x,y) (((x) < (y)) ? (x) : (y))

#define max(x,y) (((x) > (y)) ? (x) : (y))


#endif