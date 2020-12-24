#include "knn_v0.h"
#include "auxlib.h"


knnresult kNN(double * X, double * Y, int n, int m, int d, int k)
{


    /***************************************
     * Set up Environment
     ***************************************/

    // // If V0 is in V1 mode, overwrite settings
    // if(_MODE_V1_RUNNING)
    // {
    //     _TIMER_PRINT = false;
    // }

    /***************************************/


    // Start Timer
    time_t t;
    srand((unsigned) time(&t));
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);


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
        // Fix minus Zeroes (???)
        if(D2T[i]<0.0) D2T[i]=0;

        D2T[i] = (double)sqrt(D2T[i]);
    }



    if(_DIST_PRINT)
    {
        printf("--- C ---\n");

        for(int i=0; i<n; i++)
        {
            for(int j=0; j<m; j++)
            {

                if(_MODE_V1_RUNNING && i==j)
                    D2T[i*n+j] = 0;
                printf("%.2f ", D2T[i*n+j] );

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


    if(!_MODE_V1_RUNNING)
        free(D2T);


    /**
     * Result
     **/
    if(_KNN_PRINT&&!_MODE_V1_RUNNING)
    {
        print_res(res);
    }



    // Stop Timer
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    float delta_us = (float) ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000)/ (1000000);
    
    if(_TIMER_PRINT&&!_MODE_V1_RUNNING)
    {
        printf(" > V0 took %f s\n", delta_us);
    }


    return res;
}

