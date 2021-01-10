// PDS Project 2
// 
// Copyright 2021 Ⓒ K. Chatzis
// kachatzis <at> ece.auth.gr

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



/*********************************************
 *       Runtime Options & Debugging
 *********************************************/

/**
 * Parallelize index offsetting function instead
 * of running sequential. Not Optimal for few k.
 **/
// #define CILK_SPAWN_INDEX_OFFSETTING


/** 
 * Uses bins for a certain number of final
 * leafs. For V2 measurements this should be
 * defined.
 **/
#define USE_LEAF_BINS


/**
 * This definition converts the code
 * to sequential but fixes the following issue.
 * 
 * TODO: The cilk_for keyword makes the whole file's 
 * documentation unreadable on 'VS Code'. 
 * 
 * Remove before Timing!!!
 */
// #define cilk_for for



/**
 * @brief Max number of leafs to add 
 * in each search bin.
 * This should be initialized according
 * to runtime parameters.
 */
int _bin_leaf_count_max ;

// Maximum Bin Size
#define BIN_LEAF_COUNT_MAX 210
// Times without other improvements:
// N: 10^4, D: 12, No Cilk_For, No spawning
// -75 - 1.83
// 150 - 1.80
// 210 - 1.76 
// 300 - 1.96



knnresult distrAllkNNVPT(double * X, int n_all, int d, int k);



/**
 * Global Variables
 **/

int *idx_arr;       // VP index array
double *dist_arr;   // VP distance array
double *Y;          // Working Query
int N;              // Local Corpus Point count
int D;              // Vector Dimensions
int K;


/**
 * @brief VP Tree Struct. 
 * 
 *  @param vp Vantage point's Coordinates
 *  @param md Vantage point's median distance
 *  @param idx  Vantage point index
 *  @param inner Left/Inner VP Tree
 *  @param outer Right/Outer VP Tree
 *  @param bin Leaf Bin
 *  @param bin_size Size of Leaf Bin
 *  @param bin_idx Indices of bin
 * 
 *  Note: bin is initialized only for leaf nodes. 
 **/
typedef struct VPTree VPTree;
struct VPTree
{
    double *vp;
    double md; 
    int idx;   

    VPTree *inner; 
    VPTree *outer;

    double *bin;
    int *bin_idx;
    int bin_size;
};


/**
 *  Query point struct
 *  @param idx Index of point
 *  @param coord Array of vectors for each dimension (D)
 *  @param ndist Array of K nearest neighbor distances
 *  @param nidx Array of K nearest neighbor indexes
 *  @param best Closest neighbor's distance
 *  @param worst Furthest neighbor's distance
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



/**
 *  Return inner VP Tree
 *  @param T A VP Tree
 *  \return The vantage-point subtree
*/
VPTree * inner_T(VPTree * T);


/**
 *  Return outer VP Tree
 *  @param T A VP Tree
 *  \return The vantage point subtree
*/
VPTree * outer_T(VPTree * T);


/**
 *  Return median of VP Tree
 *  @param T A VP Tree
 *  \return The median distance
*/
double md_T(VPTree * T);


/**
 *  Return VP coordinates of a Tree
 *  @param T A VP Tree
 *  \return The coordinates [d-dimensional vector]
*/
double * vp_T(VPTree * T);


/**
*  Return VP index
*   @param T A VP Tree
*   \return The index to the input vector of data points
*/
int idx_T(VPTree * T);

/**
 * Checks if a VP Tree has children
 *  @param T A VP Tree
 *  \return True if the VP Tree is a leaf
*/
bool is_leaf(VPTree* T); 


/**
 *  Deletes VP Tree and all of it's children.
 *  @param T A VP Tree
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
void search_work_T(VPTree* T, Point* p);


/*
    Searches the VP Tree for the k nearest 
    neighbors.
    @param T A VP Tree
    @param p The Query point
*/
void search_vpt(VPTree* T, Point* p);


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


/**
 * @brief Prints the structure of a VP Tree.
 * 
 * @param T VP Tree
 * @param depth Used for recursive printing, 
 * should be set to zero.
 * 
 * @return nothing of use
 */
void print_vpt(VPTree* T, int depth);


/**
 * @brief Calculates Square Power of a double
 * 
 * @param x a double
 * @return double 
 */
inline double sqrd(double x) {return x*x;}


#endif  //#ifndef knn_v2_h_