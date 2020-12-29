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



void print_res(knnresult knn)
{
        for(int i=0; i<knn.m; i++)
        {
            printf("%d\nindx:", i);
            for(int j=0; j<knn.k; j++)
            {
                printf("%6d,", knn.nidx[i*knn.k+j]);
            }
            printf("\ndist:");
            for(int j=0; j<knn.k; j++)
            {
                printf("%6.2f,", (double) (knn.ndist[i*knn.k+j]));
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


void collect_n_propagate_knn(knnresult* res, knnresult* fin, int node_id, int receive_node, int send_node, int result_count)
{
    printf("in (%d) res[0].nidx[0] = %d\n",tmp_node_id ,res[1].nidx[2]);

    return;
}