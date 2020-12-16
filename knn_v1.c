#include "knn_v1.h"


#define _DIST_PRINT_VAR "DIST_PRINT"
#define _KNN_PRINT_VAR "KNN_PRINT"
#define _TIMER_PRINT_VAR "TIMER_PRINT"

// MPI Default Buffer Size
#define _MPI_COM_COUNT 1


// MPI Buffer Indexes
#define MPI_COM_STATUS_INDX 0    // Status
#define MPI_COM_DATA_MODE_INDX 1 // Data Mode
#define MPI_COM_DATA_INDX 2      // Data


// MPI Statuses
#define MPI_STATUS_INITIATING 0 // Initiating Matrices
#define MPI_STATUS_MEASURING 1  // Gathering
#define MPI_STATUS_SEARCHING 2  
#define MPI_STATUS_GATHERING 3



knnresult distrAllkNN(double * X, double * Y, int n, int m, int d, int k)
{


    /***************************************
     * Set up Environment Parameters
     * and debugging modes
     ***************************************/

    bool _knn_print = false;
    bool _dist_print = false;
    bool _timer_print = false;


    char *s_tmp;

    s_tmp = getenv( _KNN_PRINT_VAR );
    _knn_print = (s_tmp!=NULL)? ( strchr(s_tmp,'1')!=NULL? true : false  ) : false;
    // free(s_tmp);

    s_tmp = getenv( _DIST_PRINT_VAR );
    _dist_print = (s_tmp!=NULL)? ( strchr(s_tmp,'1')!=NULL? true : false  ) : false;
    // free(s_tmp);

    s_tmp = getenv( _TIMER_PRINT_VAR );
    _timer_print = (s_tmp!=NULL)? ( strchr(s_tmp,'1')!=NULL? true : false  ) : false;

    /*************************************/


    // Start Timer
    time_t t;
    srand((unsigned) time(&t));
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);


    // Initiate MPI
    MPI_Init(NULL, NULL);

    // Get Cluster Size
    int cluster_size;
    MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);

    // Get Process ID
    int mpi_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    

    




    // MPI_Recv(void *buff,
    //                int count,
    //                MPI_Datatype type,
    //                int source,
    //                int tag,
    //                int comm,
    //                MPI_Request *request )

    // MPI_Isend(void *buff,
    //                int count,
    //                MPI_Datatype type,
    //                int dest,
    //                int tag,
    //                int comm,
    //                MPI_Request *request )

    


    // double * XT;    // X Transpose
    // mat_transpose(&X, &XT, n, d);

    // double * YT;    // Y Transpose
    // mat_transpose(&Y, &YT, m, d);


    /**************************
     *       X .* X
     *        (NxD)
     **************************/
    double * X_X = (double *) malloc(n * d *sizeof(double));
    if(X_X==NULL)
    {
        printf("Failed allocating memory [X_X] (knn_v0).\n");
        exit(EXIT_FAILURE);
    }

    cilk_for(int i=0; i<n*d; i++)
    {
        X_X[i] = X[i] * X[i];
    }
    


    /***************************************
     *           D1   (Nx1)
     * for i=1:n
     *   temp = 0;
     *     for j=1:d
     *       temp = temp + XX(i,j);
     *   D1(i,1) = temp; 
     **************************************/
    double * D1 = malloc(n*sizeof(double));
    if(D1==NULL)
    {
        printf("Failed allocating memory [D1] (knn_v0).\n");
        exit(EXIT_FAILURE);
    }
    cilk_for(int i=0; i<n; i++)
    {
        double d_tmp = 0.0;
        for(int j=0; j<d; j++)
        {
            d_tmp += (double) X_X[d*i+j];  // tmp += X(i,j)
        }
        D1[i] = d_tmp;
    }

    // X_X is not needed any more
    free(X_X);


    // printf("--- D1 ---\n");
    // for(int i=0; i<n; i++){
    //     printf("%d ", (int)D1[i]);
    // }
    // printf("\n\n");



    /******************************
     *       D2 = X * Y.'   
     *          (NxM)
     ******************************/

    double * D2 = malloc(m*n*sizeof(double));

    
    /**
     * dgemm:  C = α . A x B + β . C
     * Transpose: CblasNoTrans / CblasTrans
     **/
    // NxD * DxM
    cblas_dgemm(
        CblasRowMajor,  // Store Order
        CblasNoTrans,   // A Transpose
        CblasTrans,     // B Transpose
        n, m, d,        // M, N, K of MxK * KxN
        -2.0,           // alpha
        X,              // A matrix
        d,              // K
        Y,              // B matrix
        d,              // N
        0,              // beta
        D2,             // C matrix
        m               // N
    );


    // printf("--- D2 ---\n");

    // for(int i=0; i<n; i++)
    // {
    //     // printf("%d: ( ", i);
    //     for(int j=0; j<m; j++)
    //     {

    //         printf("%d ", (int)mat_read_ij(&D2, i, j, m) );

    //     }

    //     printf("; \n");
        
    // }

    // printf("\n");




    /***************************
     *     D2 <-- D12  (NxM)
     *      (Nx1) + (NxM)
     ***************************/

    cilk_for(int i=0; i<n; i++)
    {
        double d_tmp = D1[i];
        cilk_for(int j=0; j<m; j++)
        {
            D2[i*m +j] += d_tmp;    // D2(i,j) += D1(i,1)
        }
    }



    // D1 can be freed here
    free(D1);

    // D12 is ready (in D2)


    /**************************
     *       (Y .* Y).'
     *        (DxM)
     **************************/

    double * Y_Y = malloc(m*d* sizeof(double));
    if(Y_Y==NULL)
    {
        printf("Failed allocating memory [Y_Y] (knn_v0).\n");
        exit(EXIT_FAILURE);
    }
    cilk_for(int i=0; i<m*d; i++)
    {
        Y_Y[i] = Y[i] * Y[i];
    }


    // double * Y_Y_T = malloc(d*m*sizeof(double));
    // if(Y_Y_T==NULL)
    // {
    //     printf("Failed allocating memory [Y_Y_T] (knn_v0).\n");
    //     exit(EXIT_FAILURE);
    // }

    // mat_transpose(&Y_Y, &Y_Y_T, m, d);

    // Y_Y from D3 is not needed
    // free(Y_Y);


    /***************************************
     *           D3   (1xM)
     *  for i=1:m
     *    temp = 0;
     *    for j=1:d
     *      temp = temp + YY(j,i);
     *    D3(1,i) = temp; 
     **************************************/

    double * D3 = malloc(m* sizeof(double));
    if(D3==NULL)
    {
        printf("Failed allocating memory [D3] (knn_v0).\n");
        exit(EXIT_FAILURE);
    }
    
    /**
     * Instead of doing an addition for
     * each row, we will be doing the same 
     * for each column. (~10% Improvement)
     **/
    cilk_for(int i=0; i<m; i++)
    {
        double d_tmp = 0.0;
        for(int j=0; j<d; j++)
        {
            d_tmp += (double) Y_Y[i*d+j];
        }
        D3[i] = d_tmp;
    }


    // printf("--- D3 ---\n");
    // for(int i=0; i<m; i++){
    //     printf("%d ", (int)D3[i]);
    // }
    // printf("\n\n");

    free(Y_Y);


    /***************************
     *     D2 <-- D23  (NxM)
     *      (1xM) + (NxM)
     ***************************/
    double * D2T = (double *) malloc(n*m*sizeof(double));
    mat_transpose(&D2, &D2T, n, m);
    
    free(D2);

    /** 
     * Instead of summing app for every
     * row on the Y_Y Transposed, we
     * will be adding to each column
     * of the original Y_Y.
     **/
    cilk_for(int i=0; i<m; i++)
    {
        double d_tmp = D3[i];
        cilk_for(int j=0; j<n; j++)
        {
            D2T[i*n + j] += d_tmp;
        }
    }


    // Clean Up
    free(D3);



    /*************************
     *     sqrt(D2)   (NxM)
     *************************/
    cilk_for(int i=0; i<n*m; i++)
    {
        D2T[i] = (double)sqrt(D2T[i]);
    }



    if(_dist_print)
    {
        printf("--- C ---\n");

        for(int i=0; i<n; i++)
        {
            for(int j=0; j<m; j++)
            {

                printf("%f ", (double) mat_read_ij(&D2T, j, i, n) );

            }

            printf("; \n");
        
        }

        printf("\n");
    }


    /**
     * C is ready
     **/

	knnresult res;
	res.k = k;
	res.m = m;
	res.nidx  = (int *)   malloc(m*k*sizeof(int));
	res.ndist = (double *)malloc(m*k*sizeof(double));
	
	// double * CT = malloc(n*m*sizeof(double));
    // if(CT==NULL)
    // {
    //     printf("Failed allocating memory [CT] (knn_v0).\n");
    //     exit(EXIT_FAILURE);
    // }

	// mat_transpose(&D2, &CT, n, m);
    // free(D2);


    /**
     * For each query point, spawn a 
     * thread and execute the knn search.
     * The results are returned directly 
     * into the correct places within the 
     * resulting matrices.
     **/
	cilk_for(int mm=0; mm<m; mm++)
	{
		aux_sort_idx(&D2T, &res.nidx, &res.ndist, n, m, mm, k);
	}
    

    // free(CT);



    /**
     * Result
     **/
    if(_knn_print==1)
    {
        printf("\n--- RES ---");
        for(int i=0; i<m; i++)
        {
            printf("\nindx: ", m);
            for(int j=0; j<k; j++)
            {
                printf("%d, ", (int) *(res.nidx+i*k+j));
            }
            printf("\ndist: ");
            for(int j=0; j<k; j++)
            {
                printf("%f, ", (double) *(res.ndist+i*k+j));
            }

        }
    }



    // Stop Timer
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    float delta_us = (float) ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000)/ (1000000);
    
    if(_timer_print)
    {
        printf(" > V0 took %f s\n", delta_us);
    }




    // Finalize MPI
    MPI_Finalize();

    return res;
}





// /** 
//  * Auxiliary
//  **/
// int aux_sort_idx(double** C, int** nidx, double** ndist, int N, int M, int m, int k)
// {
// 	// M: all query points
// 	// m: current query point
// 	// N: all corpus points
// 	// k: neighbors
// 	// C: distances in m*n format (Transposed)
// 	// nidx: indexes m*k
// 	// ndist: distances m*k
	
// 	// Fill the distances and indexes
// 	double * distm 	= (double *) malloc(N*sizeof(double));
// 	int * idxm 		= (int *)    malloc(N*sizeof(int));

// 	for(int i=0; i<N; i++)
// 	{
// 		distm[i] = (double) *((*C)+ N * m + i);
// 		idxm[i] = (int) i;
// 	}
	
// 	aux_mergeSort(&distm, &idxm, N);

//     for(int i=0; i<k; i++)
//     {
//         *( *ndist + k*m + i ) = (double) distm[i];
//         *( *nidx  + k*m + i ) = (int)     idxm[i];
//     }

// 	return 0;
	
// }







// /*************
//  *  Sorting 
//  *************/


// // in-place merge sort
// void aux_mergeSort(double ** _I, int ** _J, int len)
// {
//     if (len <= 1)
//     {
//         return;
//     }


//     double *left_I = (double *) cilk_spawn aux_slice_d(_I, 0, len / 2 + 1);
//     double *right_I = (double *)cilk_spawn aux_slice_d(_I, len / 2, len);

//     int *left_J = (int *) cilk_spawn aux_slice_i(_J, 0, len / 2 + 1);
//     int *right_J = (int *) cilk_spawn aux_slice_i(_J, len / 2, len);

//     cilk_sync;


//     cilk_spawn aux_mergeSort(&left_I, &left_J, len / 2);
//     cilk_spawn aux_mergeSort(&right_I, &right_J, len - (len / 2));

//     cilk_sync;

//     aux_merge(_I, _J, &left_I, &left_J, &right_I, &right_J, len / 2, len - (len / 2));

// }

// int * aux_slice_i(int **arr, int start, int end)
// {
//     int *result = (int *) malloc((end - start) * sizeof(int));
//     if(result==NULL)
//     {
//         printf("Failed Allocating memory (aux_slice_i).\n");
//     }
//     int i;
//     for (i = start; i < end; i++)
//     {
//         result[i - start] = (int) *(*arr +i);
//     }
//     return result;
// }

// double * aux_slice_d(double **arr, int start, int end)
// {
//     double *result = (double *) malloc((end - start) * sizeof(double));
//     if(result==NULL)
//     {
//         printf("Failed Allocating memory (aux_slice_d).\n");
//     }
//     int i;
//     for (i = start; i < end; i++)
//     {
//         result[i - start] = (double) *(*arr +i);
//     }
//     return result;
// }


// void aux_merge( double ** _I, int ** _J, 
//     double ** left_I, int ** left_J, 
//     double ** right_I, int ** right_J, 
//     int leftLen, int rightLen)
// {

//     int i, j;

//     i = 0;
//     j = 0;
//     while(i < leftLen && j < rightLen)
//     {
//         if ( *(*left_I +i) < *(*right_I +j) ) 
        
//         {
//             *(*_I + i + j) = *(*left_I +i);
//             *(*_J + i + j) = *(*left_J+ i);
//             i++;
//         }
//         else
//         {
//             *(*_I +i + j) = *(*right_I + j);
//             *(*_J+ i + j) = *(*right_J + j);
//             j++;
//         }

//     }

//     for(; i < leftLen; i++)
//     {
//         *(*_I +i + j) = *(*left_I + i);
//         *(*_J +i + j) = *(*left_J + i);
//     }
//     for(; j < rightLen; j++)
//     {
// 		*(*_I +i + j) = *(*right_I +j);
//         *(*_J +i + j) = *(*right_J +j);
//     }

//     free(*left_I);
//     free(*right_I);
//     free(*left_J);
//     free(*right_J);

// }



