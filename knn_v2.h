#ifndef knn_v2_h_
#define knn_v2_h_

#include "auxlib.h"
#include "mpi_wrapper.h"



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



knnresult distrAllkNNVPT(double * X, int n_all, int d, int k);



/**
 * Global Variables
 **/

int *idArr;     // VP index array
double *distArr;// VP distance array
double *Y;      // Working Query
int N;          // Local Corpus Point count
int D;          // Vector Dimensions
int K;


inline double sqrd(double x) {return x*x;}


/*!
    VP Tree Struct
    @param vp Vantage point's Coordinates
    @param md Vantage point's median distance
    @param idx  Vantage point index
    @param inner Left/Inner VP Tree
    @param outer Right/Outer VP Tree

*/
typedef struct VPTree VPTree;
struct VPTree
{
    double *vp; //the vantage point
    double md;  //the median distance of the vantage point to the others
    int idx;    //the index of the vantage point in the original set
    VPTree *inner;
    VPTree *outer;
};


/*!
    Query point struct
    @param idx Index of point
    @param coord Array of vectors for each dimension (D)
    @param ndist Array of K nearest neighbor distances
    @param nidx Array of K nearest neighbor indexes
    @param best Closest neighbor's distance
    @param worst Furthest neighbor's distance
*/
typedef struct Point Point;
struct Point
{
    int idx;
    double *coord;
    double *ndist;
    int *nidx;
    // double best;
    // double worst;
};



/*!
    Return inner VP Tree
    @param T A VP Tree
    \return The vantage-point subtree
*/
VPTree * inner_T(VPTree * T);


/*!
    Return outer VP Tree
    @param T A VP Tree
    \return The vantage point subtree
*/
VPTree * outer_T(VPTree * T);


/*!
    Return median of VP Tree
    @param T A VP Tree
    \return The median distance
*/
double md_T(VPTree * T);


/*!
    Return VP coordinates of a Tree
    @param T A VP Tree
    \return The coordinates [d-dimensional vector]
*/
double * vp_T(VPTree * T);


/*!
    Return VP index
    @param T A VP Tree
    \return The index to the input vector of data points
*/
int idx_T(VPTree * T);

/*!
    Checks if a VP Tree has children
    @param T A VP Tree
    \return True if the VP Tree is a leaf
*/
bool is_leaf(VPTree* T); 


/*!
    Deletes VP Tree and all of it's children.
    @param T A VP Tree
*/
void delete_T(VPTree* T);


/**
 * @brief Adds the Vantage point in the list of neighbors,
 * if the distance is suitable. The neighbors are always 
 * sorted on ascending distance.
 * 
 * @param T VP Tree
 * @param p Point being queried
 */
void search_work(VPTree* T, Point* p);


/*
    Searches the VP Tree for the k nearest 
    neighbors.
    @param T A VP Tree
    @param p The Query point
*/
void searchVPT(VPTree* T, Point* p);

/**
 * @brief Finds distance between two vectors. The
 * dimensions of the vectors are taken through
 * the global variable D.
 * 
 * @param p1 First Vector
 * @param p2 Second Vector
 * @return double 
 */
double distance(double* p1, double* p2);


#endif  //#ifndef knn_v2_h_