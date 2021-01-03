#include "auxlib.h"


/**
 * Each external variable should
 * be defined both in .c and .h files.
 **/

/**
 * True if V1 is running.
 * Default: False.
 **/
extern bool _MODE_V1_RUNNING = false ;

/**
 * True if V2 is running.
 * Default: False.
 **/
extern bool _MODE_V2_RUNNING = false ;

/**
 * Print kNN Neighbors result.
 * Default: False.
 **/
extern bool _KNN_PRINT = false ;

/**
 * Print calculated distance matrix.
 * Default: False.
 **/
extern bool _DIST_PRINT = false ;

/**
 * Print timing information.
 * Default: true.
 **/
extern bool _TIMER_PRINT = true ;



/** 
 * Auxiliary
 **/
int aux_sort_idx(double** C, int** nidx, double** ndist, int N, int M, int m, int k)
{
	// M: all query points
	// m: current query point
	// N: all corpus points
	// k: neighbors
	// C: distances in m*n format (Transposed)
	// nidx: indexes m*k
	// ndist: distances m*k
	
	// Fill the distances and indexes
	double * distm 	= (double *) malloc(N*sizeof(double));
	int * idxm 		= (int *)    malloc(N*sizeof(int));

	for(int i=0; i<N; i++)
	{
		distm[i] = (double) *((*C)+ N * m + i);
		idxm[i] = (int) i;
	}
	
	aux_mergeSort(&distm, &idxm, N);

    for(int i=0; i<k; i++)
    {
        *( *ndist + k*m + i ) = (double) distm[i];
        *( *nidx  + k*m + i ) = (int)     idxm[i];
    }

	return 0;
	
}




/*************
 *  Sorting 
 *************/


// in-place merge sort
void aux_mergeSort(double ** _I, int ** _J, int len)
{
    if (len <= 1)
    {
        return;
    }


    double *left_I = (double *) cilk_spawn aux_slice_d(_I, 0, len / 2 + 1);
    double *right_I = (double *)cilk_spawn aux_slice_d(_I, len / 2, len);

    int *left_J = (int *) cilk_spawn aux_slice_i(_J, 0, len / 2 + 1);
    int *right_J = (int *) cilk_spawn aux_slice_i(_J, len / 2, len);

    cilk_sync;


    cilk_spawn aux_mergeSort(&left_I, &left_J, len / 2);
    cilk_spawn aux_mergeSort(&right_I, &right_J, len - (len / 2));

    cilk_sync;

    aux_merge(_I, _J, &left_I, &left_J, &right_I, &right_J, len / 2, len - (len / 2));

}

int * aux_slice_i(int **arr, int start, int end)
{
    int *result = (int *) malloc((end - start) * sizeof(int));
    if(result==NULL)
    {
        printf("Failed Allocating memory (aux_slice_i).\n");
    }
    int i;
    for (i = start; i < end; i++)
    {
        result[i - start] = (int) *(*arr +i);
    }
    return result;
}

double * aux_slice_d(double **arr, int start, int end)
{
    double *result = (double *) malloc((end - start) * sizeof(double));
    if(result==NULL)
    {
        printf("Failed Allocating memory (aux_slice_d).\n");
    }
    int i;
    for (i = start; i < end; i++)
    {
        result[i - start] = (double) *(*arr +i);
    }
    return result;
}


void aux_merge( double ** _I, int ** _J, 
    double ** left_I, int ** left_J, 
    double ** right_I, int ** right_J, 
    int leftLen, int rightLen)
{

    int i, j;

    i = 0;
    j = 0;
    while(i < leftLen && j < rightLen)
    {
        if ( *(*left_I +i) < *(*right_I +j) ) 
        
        {
            *(*_I + i + j) = *(*left_I +i);
            *(*_J + i + j) = *(*left_J+ i);
            i++;
        }
        else
        {
            *(*_I +i + j) = *(*right_I + j);
            *(*_J+ i + j) = *(*right_J + j);
            j++;
        }

    }

    for(; i < leftLen; i++)
    {
        *(*_I +i + j) = *(*left_I + i);
        *(*_J +i + j) = *(*left_J + i);
    }
    for(; j < rightLen; j++)
    {
		*(*_I +i + j) = *(*right_I +j);
        *(*_J +i + j) = *(*right_J +j);
    }

    free(*left_I);
    free(*right_I);
    free(*left_J);
    free(*right_J);

}



void print_res(knnresult* knn)
{

    for(int i=0; i<knn->m; i++)
    {
        printf("%d\nindx:", i);
        for(int j=0; j<knn->k; j++)
        {
            printf("%6d,", knn->nidx[i*(knn->k)+j]);
        }
        printf("\ndist:");
        for(int j=0; j<knn->k; j++)
        {
            printf("%6.2f,", (double) (knn->ndist[i*(knn->k)+j]));
        }

        printf("\n");
    }

}



_runtime startup(int argc, char** argv)
{


    char _x_filename[1024];

    bool _x_value_included = true;
    bool _x_transpose = false;

    char _y_filename[1024];

    bool _y_value_included = true;
    bool _y_transpose = false;

    int _n = 1000;
    int _m = 100;
    int _d = 10;
    int _k = 5;

    bool _rand = false;



    /***************************************
     * Set up Environment Parameters
     * and debugging modes
     ***************************************/

    char *s_tmp;

    s_tmp = getenv( _KNN_PRINT_VAR );
    _KNN_PRINT = (s_tmp!=NULL)? ( strchr(s_tmp,'1')!=NULL? true : false  ) : false;
    // free(s_tmp);

    s_tmp = getenv( _DIST_PRINT_VAR );
    _DIST_PRINT = (s_tmp!=NULL)? ( strchr(s_tmp,'1')!=NULL? true : false  ) : false;
    // free(s_tmp);

    s_tmp = getenv( _TIMER_PRINT_VAR );
    _TIMER_PRINT = (s_tmp!=NULL)? ( strchr(s_tmp,'1')!=NULL? true : false  ) : false;

    /***************************************/




    /** 
    * For each argument given,
    * find the respective configuration
    * parameter and change it
    **/
    for(int i=0; i<argc; i++){


        int _tmp_int = NULL;   // Temporary placeholder
        char _tmp_str[1024];


        /**
         * " -transpx -TraspX -xTransp 
         * -XTransp -xtranspose -xTranspose "
         */
        if( strcmp(argv[i],"-transpx")==0 )
        {
            _x_transpose = true;
            continue;
        }


        /**
         * " -transpy -TraspY -yTransp 
         * -YTransp -ytranspose -yTranspose "
         */
        if( strcmp(argv[i],"-transpy")==0 )
        {
            _y_transpose = true;
            continue;
        }


        /**
         * " -t_ "
         */
        if(sscanf(argv[i], "-t%d", &_tmp_int))
        {
            // if(_tmp_int>0 && _tmp_int<64)
            // {
            //     __threads = _tmp_int;
            //     continue;
            // }
            printf("-t Parameter Deprecated: Please use 'export CILK_NWORKERS=_' instead.\n");
        }



        /**
         * " -k_ "
         */
        if(sscanf(argv[i], "-k%d", &_tmp_int))
        {
            if(_tmp_int>0)
            {
                _k = _tmp_int;
                continue;
            }
        }


        /**
         * " -m_ "
         */
        if(sscanf(argv[i], "-m%d", &_tmp_int))
        {
            if(_tmp_int>0)
            {
                _m = _tmp_int;
                continue;
            }
        }


        /**
         * " -n_ "
         */
        if(sscanf(argv[i], "-n%d", &_tmp_int))
        {
            if(_tmp_int>0)
            {
                _n = _tmp_int;
                continue;
            }
        }


        /**
         * " -d_ "
         */
        if(sscanf(argv[i], "-d%d", &_tmp_int))
        {
            if(_tmp_int>0)
            {
                _d = _tmp_int;
                continue;
            }
        }


        /**
         * " -rand "
         */
        if( strcmp(argv[i],"-rand")==0 )
        {
            _rand = true;
            continue;
        }


        // /**
        //  * " -x< Path > "
        //  */
        // if(sscanf(argv[i], "-x%s", _tmp_str))
        // {
        //     strcpy(_x_filename, _tmp_str);
        //     printf("xfile 1 %s\n", _x_filename);
        //     continue;
        // }


        /**
         * " -x < Path > "
         */
        if(strcmp(argv[i],"-x")==0)
        {
            if(i<argc-1)
            {
                
                strcpy(_x_filename, argv[i+1]);
                continue;
            }
        }


        // /**
        //  * " -y< Path > "
        //  */
        // if(sscanf(argv[i], "-y%s", _tmp_str))
        // {
        //     strcpy(_y_filename, _tmp_str);
        //     printf("yfile 1 %s\n", _y_filename);
        //     continue;
        // }


        /**
         * " -y < Path > "
         */
        if(strcmp(argv[i],"-y")==0)
        {
            if(i<argc-1)
            {
                strcpy(_y_filename, argv[i+1]);
                continue;
            }
        }


    }



    /**
    * Configure threads on system runtime level
    **/

    // Cilk:
    // char _tmp_str[4];
    // sprintf(_tmp_str, "%d", __threads);
    // __cilkrts_set_param("nworkers", _tmp_str);

    // // OpenMP:   
    // omp_set_num_threads(__threads);
    


    /**
     * Corpus Points
     * n-by-d
     **/
    double * X; 

    /**
     * Query Points
     * m-by-d
     **/
    double * Y;

    /**
     * Corpus points count
     **/
    int n = (_rand) ? _n : 0;

    /**
     * Query points count
     **/
    int m = (_rand) ? _m : 0;

    /**
     * Dimension count
     * 
     * ! It has to be initialized as 0,
     * otherwise mmarket_import might keep an 
     * arbitrary number of dimensions.
     **/
    int d = (_rand) ? _d : 0;
    
    /**
     * neighbours count
     **/
    int k = _k;



    if(!_rand)
    {

        /*******************
         ** Import Matrix **
        ********************/

        mmarket_import(_x_filename, &X, &n, &d, _x_value_included, _x_transpose, true); // Import MM

        printf("n: %d, d: %d\n", n, d);

        mmarket_import(_y_filename, &Y, &m, &d, _y_value_included, _y_transpose, true); // Import MM

        printf("m: %d, d: %d\n", m, d);

    }
    else
    {
        

        /*******************
        ** Create Matrix **
        ********************/

        // Use current time as seed for random generator
        if(RAND_SEED)
            srand(time(0));

        // X
        X = (double *) malloc(n*d*sizeof(double));
        if(X==NULL)
        {
            printf("Failed Allocating memory (main).\n");
        }
        for(int i=0; i<n; i++)
        {
            for(int j=0; j<d; j++)
            {
                X[i*d+j] = (double) rand()/RAND_MAX*MAX_VECTOR;
            }
        }


        // Y
        Y = (double *) malloc(n*d*sizeof(double));
        if(Y==NULL)
        {
            printf("Failed Allocating memory (main).\n");
        }
        for(int i=0; i<m; i++)
        {
            for(int j=0; j<d; j++)
            {
                Y[i*d+j] = (double) rand()/RAND_MAX*MAX_VECTOR;
            }
        }

    }


    if(k>n)
    {
        printf("Illegal Parameters (k>n).\n");
        exit(EXIT_FAILURE);
    }

    if(d<1||d>200)
    {
        printf("Illegal Parameters (d<1 or d>200).\n");
        exit(EXIT_FAILURE);
    }

    if(n<1 || m<1)
    {
        printf("Illegal Parameters (n or m <0).\n");
        exit(EXIT_FAILURE);
    }


    if(MATRIX_PRINT)
    {
        printf("\n");

        printf("--- X ---\n");

        for(int i=0; i<n; i++)
        {
            // printf("", i);
            for(int j=0; j<d; j++)
            {

                printf("%7.2f", (double)mat_read_ij(&X, i, j, d) );

            }

            printf("; \n");
            
        }


        //     mat_transpose(&X, &Y, n, d);

        printf("\n");

        // printf("--- Y ---\n");

        // if(!V0_USE_X_AS_Y)
        // {
        //     for(int i=0; i<m; i++)
        //     {
        //         // printf("%d: ( ", i);
        //         for(int j=0; j<d; j++)
        //         {

        //             printf("%7.2f", mat_read_ij(&Y, i, j, d) );

        //         }

        //         printf("; \n");
                
        //     }
        // }else{
        //     printf(" [Same as X]\n");
        // }

    }


    _runtime r;
    r.n = n;
    r.m = m;
    r.k = k;
    r.d = d;
    r.X = X;
    r.Y = Y;

    return (_runtime) r;

}



extern int tmp_node_id=-1;


void collect_n_propagate_knn(knnresult* res, knnresult* fin, int node_id, int node_receive, int node_send, int result_count)
{
    /**
     * Data Collection:
     * 
     * Worker nodes should receive the
     * kNN result of the previous node
     * (except node=1, doesn't receive)
     * and compare/update with the local
     * results. Next, the receiving nodes
     * should send the result to the next
     * node in the ring.
     * 
     * Master node only receives and compares.
     * All worker nodes have finished
     * after they have sent their kNN
     * results.
     * 
     * The result tranfered is a single 
     * kNN result, therefore a loop is
     * used to account for all batches of
     * Query points.
     **/


    if(fin==NULL || res==NULL)
        mpi_abort_msg("fin or res were NULL (collect_n_propagate_knn)");
    
    if(fin->m == NULL || fin->m <0)
        mpi_abort_msg("Illegal or NULL m (collect_n_propagate_knn)");

    if(fin->k == NULL || fin->k < 0)
        mpi_abort_msg("Illegal or NULL k (collect_n_propagate_knn)");
    
    int n_all = fin->m;
    int k = fin->k;

    int m_all = 0;  // Used to test points imported

    // Cluster Size = Result Count = Batch Count
    int cluster_size = result_count;    

    // Request array for receive and send
    MPI_Request mpi_request[2];

    
    // mpi_finish_local();
    // printf("n: %d, k: %d\n", n_all, k);


    // wait( (float)rand()/(float)(RAND_MAX/2) );

    // printf("******\nNode: %d:\n", node_id);
    // for(int i=0; i<result_count; i++)
    // {
    //     print_res(&res[i]);
    // }





    // Nodes that receive
    if(node_id!=1)
    {

        // For each batch
        for(int i=0; i<cluster_size; i++)
        {
            // TODO: Remove m receiving,
            // We already know m is n_all/cluster_size.
            // The last one will have locally the same
            // m plus n_all%cluster_size.
            
            // Receive M
            int m_buff[1];
            mpi_receive_data_b_t(
                MPI_MODE_KNN_COLLECTION_M, 
                m_buff, 
                1, 
                node_receive, 
                MPI_INT, 
                mpi_request
            );

            m_all += m_buff[0];
            // if(DEBUG_A)
                // printf("[%d] batch: %d, m: %d\n", node_id, i, m_buff[0]);

            // if(m_buff[0] < 0)
            // {
            //     printf("M received from node %d on node %d was illegal!\n");
            //     exit(EXIT_FAILURE);
            // }

            // TODO: Maybe execute 2 concurrent receives

            int * indx_buff = (int*)malloc(m_buff[0]*k*sizeof(int));
            mpi_receive_data_b_t(
                MPI_MODE_KNN_COLLECTION_INDX, 
                indx_buff,
                m_buff[0]*k,
                node_receive,
                MPI_INT,
                mpi_request
            );

            double * dist_buff = (double*)malloc(m_buff[0]*k*sizeof(double));
            mpi_receive_data_b_t(
                MPI_MODE_KNN_COLLECTION_DIST, 
                dist_buff,
                m_buff[0]*k,
                node_receive,
                MPI_DOUBLE,
                mpi_request
            );

            knnresult* res_tmp = malloc(sizeof(knnresult));
            res_tmp->ndist = dist_buff;
            res_tmp->nidx  = indx_buff;
            res_tmp->m     = m_buff[0];
            res_tmp->k     = k;

            cilk_spawn compare_knnresult(&res[i], res_tmp);

        }

        cilk_sync;
            
    }//end: if(node_id>1)


    // nodes that send
    if(node_id!=0)
    {
        // TODO: Move sending right after calculation
        // to make it async.

        // Send results to next node
        for(int i=0; i<cluster_size; i++)
        {
            // m
            int m_buff[] = { res[i].m };
            mpi_send_data_b_t(
                MPI_MODE_KNN_COLLECTION_M, 
                m_buff,
                1,
                node_send,
                MPI_INT,
                mpi_request
            );

            // nidx
            mpi_send_data_b_t(
                MPI_MODE_KNN_COLLECTION_INDX, 
                res[i].nidx,
                res[i].m*k,
                node_send,
                MPI_INT,
                mpi_request
            );

            // ndist
            mpi_send_data_b_t(
                MPI_MODE_KNN_COLLECTION_DIST, 
                res[i].ndist,
                res[i].m*k,
                node_send,
                MPI_DOUBLE,
                mpi_request
            );

        }
    }

    // printf("m_all:%d\n", m_all);
    // if(m_all != n_all && node_id!=1)
    //     mpi_abort_msg("Received Query Point results were wrong");

    // Nodes other than master have finished
    if(node_id>0)
    {
        mpi_send_data_wait(mpi_request);
        // mpi_receive_data_wait(mpi_request);
        return;
    }


    // Only Master node continues here

    // Concetrate results on a single knnresult

    double* knn_ndist = (double*) malloc(m_all*k*sizeof(double));
    int* knn_nidx     =    (int*) malloc(m_all*k*sizeof(int))   ;

    if(knn_ndist==NULL)
        mpi_abort_msg("Malloc Failed (knn_ndist in collect_n_propagate_knn)");
    if(knn_nidx==NULL)
        mpi_abort_msg("Malloc Failed (knn_nidx in collect_n_propagate_knn)");
            
    

    for(int i=0; i<cluster_size; i++)
    {

        int lgth = (i==cluster_size-1) ? 
            (n_all/cluster_size)+n_all%cluster_size : 
            n_all/cluster_size ;

        // wait( (float)rand()/(float)(RAND_MAX) );
        // printf("Node: %d\n", node_id);
        // printf("\ni %d\n", i);
        // printf("--indxx:\n");
        // for(int j=0; j<res[i].m; j++)
        // {
        //     printf("  ");
        //     for(int p=0; p<res[i].k; p++){
        //         printf("%d, ", res[i].nidx[j*res[i].k + p]);
        //     }
        //     printf("\n");
        // }
        // printf("\n");

        // printf("--distx:\n");
        // for(int j=0; j<res[i].m; j++)
        // {
        //     printf("  ");
        //     for(int p=0; p<res[i].k; p++){
        //         printf("%.2f, ", res[i].ndist[j*res[i].k + p]);
        //     }
        //     printf("\n");

        // }
        // printf("\n");
        // printf("--lgth:%d\n", lgth);

        // Copy batch ndist
        cilk_spawn memcpy(
            &knn_ndist[ i*k*(n_all/cluster_size) ],
            &(res[i].ndist[0]),
            res[i].k * res[i].m * sizeof(double)
        );

        // Copy batch nidx
        cilk_spawn memcpy(
            &knn_nidx[ i*k*(n_all/cluster_size) ],
            &(res[i].nidx[0]),
            res[i].k * res[i].m * sizeof(int)
        );

    }

    cilk_sync; // Sync result concetration

    // knnresult knn;
    fin->m = n_all;
    fin->k = k;
    fin->ndist = knn_ndist;
    fin->nidx = knn_nidx;

    return;
}


void knnresult_check_offset(knnresult* res, const int min, const int max)
{
    int mmin = max;
    int mmax = min;

    for(int i=0; i<res->m * res->k; i++)
    {
        if(res->nidx[i]>mmax)
            mmax = res->nidx[i];

        if(res->nidx[i]<mmin)
            mmin = res->nidx[i];
    }

    if(mmin<min || mmax>max)
        {
            printf("Offset Check Failed! mmin: %d (%d)\nmmax: %d (%d)\n", mmin, min, mmax, max);
            mpi_abort_msg("Offset Check");
        }
    
}

void knnresult_check_batch(double *Y, double *X, int length, int offset)
{
    // m = n * d
    for(int i=0; i<length; i++)
    {
        if(Y[i] != X[i+offset])
        {
            printf("Batch Check Failed! (found:%.2f, correct:%.2f , length:%d, offset:%d, i:%d)\n", Y[i], X[i+offset], length, offset, i);
            mpi_abort_msg("Batch Check");
        }
    }
}


void knnresult_offset_nidx(knnresult* res, int offset)
{
    for(int i=0; i<res->m * res->k; i++)
    {
        res->nidx[i] += offset;
    }
}


int n_per_node(int node_id, int cluster_size, int n)
{

    int ret = (n/cluster_size);
    if(node_id == cluster_size-1)
    {
        return ret + n%cluster_size;
    }
    return ret;
}


void compare_knnresult(knnresult* res, knnresult* res_)
{

    res->m = (res->m > res_->m)? res->m : res_->m ;

    int alloc = 2*res->k;


    for(int m = 0; m<res->m; m++)
    {

        double * distm 	= (double *) malloc(alloc*sizeof(double));
        int * idxm 		= (int *)    malloc(alloc*sizeof(int));

        for(int i=0; i<res->k; i++)
        {
            distm[i] = (double) res->ndist[m*res->k+i];
            idxm[i]  =    (int) res->nidx[m*res->k+i] ;
        }

        for(int i=0; i<res_->k; i++)
        {
            distm[i+res->k] = (double) res_->ndist[m*res->k+i];
            idxm[i+res->k]  =    (int) res_->nidx[m*res->k+i] ;
        }

        aux_mergeSort(&distm, &idxm, alloc);

        for(int i=0; i<res->k; i++)
        {
            res->ndist[m*res->k+i] = (double) distm[i];
            res->nidx[m*res->k+i]  =    (int) idxm[i] ;
        }

    }

}