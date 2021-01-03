#include "knn_v2.h"

#define DEBUG_VPSEARCH false
#define DEBUG_VPTREE false


/**
 * TODO: The cilk_for keyword makes the whole file's 
 * documentation unuseable. This definition converts
 * the code to sequential but fixes the issue.
 * Remove before Release!!!
 */
// #define cilk_for for



// int tmp_batch_id = -1;

void print_vpt(VPTree* T, int depth)
{
    if(T==NULL)
        return;
        
    char* pad = (char*)malloc( (depth+1)*sizeof(char) );
    int i;
    for(i=0; i<depth; i++)
    {
        pad[i] = ' ';
    }
    pad[i] = '\0';
    printf("%s%2d (%5.2f)       [%5.2f]\n", pad, T->idx,T->vp[0], T->md);
    // printf("i");
    print_vpt((T->inner), depth+2);
    // printf("o");
    print_vpt((T->outer), depth+2);
}


knnresult distrAllkNNVPT(double * X, int n_all, int d, int k)
{

    // Start Timer
    time_t t;
    srand((unsigned) time(&t));
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);


    /**
     * Set V0 to V2 mode. 
     * (Not used for anything yet?)
     **/
    _MODE_V2_RUNNING = true;

    // Initialize MPI
    int node_id, cluster_size;
    mpi_initialize(&node_id, &cluster_size);

    // Request array for receive and send
    MPI_Request mpi_request[2];


    if(cluster_size<2)
    {
        printf("There is no cluster!\n");
        mpi_abort();
    }

    if(k > n_all/cluster_size)
    {
        printf("Illegal K!\n");
        mpi_abort();
    }


    // Who I receive from 
    int node_receive = (node_id==0)? cluster_size-1 : node_id-1;

    // Who I send to
    int node_send = (node_id==cluster_size-1)? node_send = 0 : node_id+1;

    if(node_receive>=cluster_size||node_receive<0||node_send>=cluster_size||node_send<0)
        mpi_abort_msg("Illegal Node Send/Receive");


    int n = n_per_node(node_id, cluster_size, n_all);

    N = n;  // Local Corpus Points (Global)
    D = d;  // Dimensions (Global)
    K = k;  // Nearest Neighbors


    if(node_id>0)
    {
        mpi_receive_data_b(MPI_MODE_DATA_DISTRIBUTION, X, n_all*d, node_receive, mpi_request);
    }

    // All nodes, except the last one, 
    // should send the matrix to the
    // next one.
    if(node_id != cluster_size-1)
    {
        mpi_send_data_b(MPI_MODE_DATA_DISTRIBUTION, X, n_all*d, node_send, mpi_request);
    }



    double *XX = NULL;

    XX = malloc(n_all*d*sizeof(double));
    if(XX==NULL) mpi_abort_msg("Malloc Failed (XX)");
    memcpy(&XX[0], &X[0], n_all*d*sizeof(double));


    int node_offset  = node_id*(n_all/cluster_size);


    // Hold on to local data only
    memcpy(&X[0], &X[ d*node_offset ], d*n_per_node(node_id, cluster_size, n_all) *sizeof(double));
    X = realloc(X, d*n_per_node(node_id, cluster_size, n_all)*sizeof(double));
    if(X==NULL) mpi_abort_msg("Realloc Failed (X)");


    // Allocate
    // Allocation enough for the largest batch
    int yz_alloc_size = d*n_per_node(cluster_size-1, cluster_size, n_all);

    Y = (double*) malloc(yz_alloc_size*sizeof(double));
    if(Y==NULL) mpi_abort_msg("Malloc Failed (Y)");

    // Copy X to Y
    memcpy(&Y[0], &X[0], d*n_per_node(node_id, cluster_size, n_all) *sizeof(double));
    free(X);

    double* Z = (double*) malloc(yz_alloc_size*sizeof(double));
    if(Z==NULL) mpi_abort_msg("Malloc Failed (Z)");

    // wait( (float)rand()/(float)(RAND_MAX) );
    // printf("(%d) X: %f\n\n", node_id, X[0]);
    

    /**
     * TODO: Use same root VP in each node.
     **/


    /** 
     * Build VP Tree for local corpus
     **/

    VPTree *T   = (VPTree*) malloc( sizeof(VPTree)  );
    idArr       = (int*)    malloc( n_per_node(node_id, cluster_size, n_all)*sizeof(int)   );
    distArr     = (double*) malloc( n_per_node(node_id, cluster_size, n_all)*sizeof(double));

    if(T==NULL)      mpi_abort_msg("Malloc Failed (T)");
    if(idArr==NULL)  mpi_abort_msg("Malloc Failed (idArr)");
    if(distArr==NULL)mpi_abort_msg("Malloc Failed (distArr)");
    
    for(int i=0; i<n_per_node(node_id, cluster_size, n_all); i++) 
    {
        idArr[i] = i;
    }

    recursiveBuildTree(T, 0, n_per_node(node_id, cluster_size, n_all)-1);
    

    // free(idArr);
    // idArr = NULL;

    // free(distArr);
    // distArr = NULL;


    wait( (float)rand()/(float)(RAND_MAX) );
    if(DEBUG_VPTREE)
        print_vpt(T, 0);    // Print VP Tree for Debugging
        

    // kNN Results per batch
    knnresult* result = (knnresult*)malloc(cluster_size*sizeof(knnresult));
    

    // Working Batch ID
    int working_batch = node_id;
    
    // Query Point Count
    int m = n_all/cluster_size;
    m += (working_batch==cluster_size-1)? n_all%cluster_size: 0;

    // Query Point batch offset for indices
    int batch_offset;
    // wait( (float)rand()/(float)(RAND_MAX/3) );

    int zeroes_found = 0;
    for(int w=0; w<cluster_size; w++)
    {

        
        // printf("Batch %d, Node %d\n", working_batch, node_id);

        tmp_batch_id = working_batch;
        tmp_node_id = node_id;

        batch_offset = working_batch*(n_all/cluster_size);

        if(working_batch!=node_id)
        {
            mpi_send_data_wait(mpi_request);
            mpi_receive_data_wait(mpi_request);
            memcpy(&Y[0], &Z[0], m*d*sizeof(double));
        }
            

        // memcpy(&Y[0], &Z[0], m*d*sizeof(double));
        // Send Current Working Batch
        mpi_send_data_nb(MPI_MODE_QUERY_DISTRIBUTION, Y, m*d, node_send, mpi_request);

        // wte: What To Expect
        // What length to expect, depending
        // on the next batch id
        int wte = (working_batch-1<0)?cluster_size-1:working_batch-1;
        // if(wte<0)
        //     wte = cluster_size -1;
        wte = n_per_node(wte, cluster_size, n_all);

        // printf("(%d) Receiving Batch %d, Size %d\n", node_id, working_batch-1, wte);

        // Receive Next batch
        // mpi_receive_data_nb(MPI_MODE_CORPUS_DISTRIBUTION, Z, wte, node_receive, mpi_request);
        mpi_receive_data_nb(MPI_MODE_QUERY_DISTRIBUTION, Z, wte*d, node_receive, mpi_request);
        // printf("(%d) got batch %d\n", node_id, working_batch);


        

        /**
         * Set Up Result Space
         **/
        int res_int_alloc = (n_all/cluster_size + n_all%cluster_size) *k;
        result[working_batch].ndist = (double*)malloc(res_int_alloc*sizeof(double));
        result[working_batch].nidx  = (int*)malloc(res_int_alloc*sizeof(int));
        result[working_batch].m = m;
        result[working_batch].k = k;



        // For each point in the query
        cilk_for(int p=0; p<m; p++)
        {
            Point* _point = (Point*)malloc(sizeof(Point));
            _point->coord = (double*)malloc(D*sizeof(double));
            for(int dd=0; dd<D; dd++)
            {
                (_point->coord)[dd] = Y[p*D+dd];
            }
            _point->ndist = (double*)calloc(k,sizeof(double));
            _point->nidx  = (int*)malloc(k*sizeof(int));
            for(int kk=0; kk<K; kk++)
            {
                _point->nidx[kk] = -1;
            }
            _point->idx = p + batch_offset;

            wait( (float)rand()/(float)(RAND_MAX/2) );
            // if(T->vp[0] != XX[node_offset+T->idx])
            // {
                // printf("\n\n ## DIFF ## Node: %d/%d, (%5.2f: %5.2f) \n\n",node_id, working_batch,XX[node_offset+T->idx], T->vp[0]);
                if(DEBUG_VPTREE)
                    print_vpt(T, 0);    // Print VP Tree for Debugging
        
            // }
            searchVPT(T, _point);

            
            // if(_point->ndist[0]==0.0 && node_id!=working_batch)
            // {
            //     printf("\n*****CAUGHT*****\n");
            // }

            // Offset resulting indexes
            cilk_for(int i=0; i<k; i++)
            {
            //     if(_point->nidx[i] < 0 || _point->ndist[i]<0 ||
            //         _point->nidx[i] > n_all || _point->ndist[i] > 40 )
            //         mpi_abort_msg("--- Illegal ndist/nidx ---\n");
                _point->nidx[i] += node_offset;
                // printf("(%d/%d) Point: %d /%6.2f, kk: %d, idx: %d, dist:%6.2f\n", node_id, working_batch, _point->idx, _point->coord[0], i, _point->nidx[i], _point->ndist[i]);
            }


            // Add it into result[working_batch]
            cilk_for(int kk=0; kk<k; kk++)
            {
                (result[working_batch]).ndist[p*k+kk] = _point->ndist[kk];
                (result[working_batch]).nidx[p*k+kk] = _point->nidx[kk];
            }

            // for(int i=0; i<k; i++)
            // {
            //     if(result[working_batch].ndist[p*k+i]==0.0 && node_id!=working_batch)
            //     {

            //         printf("\n!!!(%d/%d) K: %d, Dist:%6.2f\n", node_id, working_batch, i, result[working_batch].ndist[p*k+i]);
            //         // mpi_abort_msg("Zero distance on different batches");
            //     }
            // }

            
            // Clean Up
            free(_point->coord);
            free(_point->ndist);
            free(_point->nidx);
            free(_point);
            _point = NULL;

        }


        // update working batch id
        working_batch--;
        if(working_batch<0)
            working_batch = cluster_size-1;

        // update batch's query point count
        m = n_all/cluster_size;
        m += (working_batch==cluster_size-1)? n_all%cluster_size: 0;

    }
// wait( (float)rand()/(float)(RAND_MAX/2) );
    // TODO: Not Required?
    mpi_send_data_wait(mpi_request);
    


    // mpi_finish_local();
    // Clean Up
    
    
    delete_T(T);
    T = NULL;


    tmp_node_id = node_id; // For Debugging auxlib, TODO: Remove

    knnresult* knn = (knnresult*)malloc(sizeof(knnresult));
    knn->m = n_all;
    knn->k = k;


    collect_n_propagate_knn(result, knn, node_id, node_receive, node_send, cluster_size);
    // wait( (float)rand()/(float)(RAND_MAX/2) );

    // sleep(cluster_size-node_id+2);

    /**
     * Worker Nodes should collect knnresults,
     * compare them, and propagate clockwise.
     **/
    if(node_id>0)
    {
        // Clean Up
        for(int i=0; i<cluster_size; i++)
        {
            free(result[i].ndist);
            free(result[i].nidx);
        }
        free(result);

        free(knn);
        
        mpi_finish_local(); // bye
    }

    // Only the Master node continues here

    // printf("-- END (%d) --\n", node_id);
    // return;
    // mpi_finish_local();
    

    // printf("m: %d, k: %d\n", knn->m, knn->k);

    // if(knn->ndist == NULL)
    //     printf("ndist NULL !!!]n");
    // printf("dist 1 1 : %f\n", (knn->ndist[1]) );
    // mpi_finish_local();

    // sleep(2);
    // printf("-------------------------\n");

    if(_KNN_PRINT)
    {
        print_res(knn);
    }



    // Stop Timer
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    float delta_us = (float) ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000)/ (1000000);
    if(_TIMER_PRINT)
    {
        printf(" > V1 took %f s\n", delta_us);
    }

    return *knn;

    mpi_abort_msg("Master node didn't return anything");

}



void search_work(VPTree* T, Point* p)
{

    if(DEBUG_VPSEARCH)
        printf("  Work, VP: %d (%5.2f), P: %d(%5.2f)\n", T->idx, T->vp[0], p->idx, p->coord[0]);
    
    double d_d = distance((double*)T->vp, (double*) p->coord); // VP-to-Point distance
    int d_i = T->idx;   // VP index

    // Auxiliary variables for swaps
    double r_d;
    int r_i;

    if(DEBUG_VPSEARCH)
    {
        printf("    Dist: %5.2f\n", (double)d_d);

        printf("    neighbors_pre:{");
        for(int i=0; i<K; i++)
        {
            printf("%2d[%5.2f],", p->nidx[i], p->ndist[i]);
        }
        printf("}\n");
    }

    if( d_d < p->ndist[K-1] || p->nidx[K-1]<0 )    // VP distance is within the neighbors
    {
        // Push point to neighbors
        
        /**
         * TODO: use Binary Search
         **/

        // For each neighbor
        for(int i=0; i<K; i++)
        {
            // If the distance is sorter
            if(p->ndist[i] > d_d || p->nidx[i]<0)
            {
                // Swap
                r_i = p->nidx[i];
                r_d = p->ndist[i];

                p->ndist[i] = d_d;
                p->nidx[i] = d_i;

                d_i = r_i;
                d_d = r_d;
            }
        }

    }

    if(DEBUG_VPSEARCH)
    {
        printf("    neighbors_pos:{");
        for(int i=0; i<K; i++)
        {
            printf("%2d[%5.2f],", p->nidx[i], p->ndist[i]);
        }
        printf("}\n");
    }

    return;
}


double distance(double* p1, double* p2)
{
    // printf("    Distance(%5.2f, %5.2f)\n", p1[0], p2[0]);
    double res = 0;
    for(int i=0; i<D; i++)
    {
        res += sqrd(p2[i] - p1[i]);
    }
    return (double)sqrt((double)res);
}


void searchVPT(VPTree* T, Point* p)
{
    
    // Check pointer are present
    if(T==NULL || p==NULL)
        return;

    wait( (float)rand()/(float)(RAND_MAX/2) );

    if(DEBUG_VPSEARCH)
        printf("(%d/%d) SearchVPT, T:%d, P:%d  (%5.2f)\n", tmp_node_id, tmp_batch_id, T->idx, p->idx, (p->coord)[0]);


    // Finds distance of p from vp.
    // If distance is less than the most distant neighbor,
    // push vantage point into p's neighbors
    search_work(T, p);
    

    if(is_leaf(T))   // T has no children
    {
        return ;
    }
    else    // T has children
    {
        
        // I should first where I belong
        // (left or right).
        // After returning, if there is no 
        // intersection, search the other side.

        if( distance((double*)T->vp, (double*)p->coord) < T->md )
        {
            if(DEBUG_VPSEARCH)
                printf("...Inner\n");

            searchVPT(T->inner, p);

            // TODO: ..
            // if( check_intersection )
            // {
                searchVPT(T->outer, p);
            // }
        }
        else /* if ( distance(T->vp, p) > T->md  )*/
        {
            if(DEBUG_VPSEARCH)
                printf("...Outer\n");
            searchVPT(T->outer, p);

            // if( check_intersection )
            // {
                searchVPT(T->inner, p);
            // }
        }


    }
    

}






void calc_distance(double *vp, int start, int end)
{
    // Coordinates as [i][j] instead of [i*j_max+j]
    double(*dataArr)[D] = (double(*)[D])Y;

    // Squared Distance of each node
    for (int i=start; i<=end; i++)
    {
        distArr[i] = sqrd(vp[0] - dataArr[idArr[i]][0]);
    }

    for (int i=start; i<=end; i++)
    {
        for (int j=1; j<D; j++)
        {
            distArr[i] += sqrd(vp[j] - dataArr[idArr[i]][j]);
        }
    }

}



void swap_d(double* a, double* b)
{
    double temp = *a;
    *a = *b;
    *b = temp;
}


void swap_i(int* a, int* b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}



// Copied
void quick_select(int kpos, double* distArr, int* idArr, int start, int end)
{
    int store=start;
    double pivot=distArr[end];
    for (int i=start; i<=end; i++)
        if (distArr[i] <= pivot)
        {
            swap_d(distArr+i, distArr+store);
            swap_i   (idArr+i,   idArr+store);
            store++;
        }        
    store--;
    if (store == kpos) return;
    else if (store < kpos) quick_select(kpos, distArr, idArr, store+1, end);
    else quick_select(kpos, distArr, idArr, start, store-1);
}



void recursiveBuildTree(VPTree* node, int start, int end)
{

    /**
     * TODO: Stop at a constant depth
     **/

    /**
     * TODO: Maye remove Square root to make it faster
     **/
    
    // Get Node ID
    node->idx = idArr[end];
    
    // Get Node Coordinates from Y
    node->vp = (double*)malloc(D*sizeof(double));
    memcpy((node->vp), &Y[node->idx*D], D*sizeof(double));
    
    if (start==end)
    {
        // Initialize
        node->inner = NULL;
        node->outer = NULL;
        node->md = 0.0; 
        return;
    }

    end--; // VP was at the end, we remove it

    // Calculate distances to all other nodes
    calc_distance(node->vp, start, end);
    
    quick_select( (start+end)/2, distArr, idArr, start, end );

    // idArr is half inner and half outer points (from
    // the perspective of the median)

    // Median is the point found at half array length
    node->md    = sqrt(distArr[ (start+end)/2 ]);

    // Inner will continue constructing from the median backwards to the start
    node->inner = (VPTree*)malloc( sizeof(VPTree) );

    // Outer will construct from the end backwards, up to the median
    node->outer = (VPTree*)malloc( sizeof(VPTree) );

    recursiveBuildTree(node->inner, start, (start+end)/2);
    if (end>start)
    {
        recursiveBuildTree(node->outer, (start+end)/2 +1, end);
    }
    else 
    {
        node->outer = NULL;
    }

    return;
}


VPTree* inner_T(VPTree* T) { return T->inner; }

VPTree* outer_T(VPTree* T) { return T->outer; }

double md_T(VPTree* T)     { return T->md; } 

double* vp_T(VPTree* T)    { return T->vp; }

int idx_T(VPTree* T)       { return T->idx; }


bool is_leaf(VPTree* T) 
{
    return (inner_T(T)==NULL&&inner_T(T)==NULL);
}



void delete_T(VPTree* T)
{

    if(T==NULL)
        return;

    cilk_spawn delete_T(T->inner);
    cilk_spawn delete_T(T->outer);

    cilk_spawn free(T);

    cilk_sync;

    return;
}



int rand_i(int min, int max)
{
    rand() % (max + 1 - min) + min;
}



