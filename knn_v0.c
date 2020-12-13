#include "knn_v0.h"



// void main()
// {
//   int i=0;
//   double A[6] = {1.0,2.0,1.0,-3.0,4.0,-1.0};
//   double B[6] = {1.0,2.0,1.0,-3.0,4.0,-1.0};
//   double C[9] = {.5,.5,.5,.5,.5,.5,.5,.5,.5};
//   cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,3,3,2,1,A, 3, B, 3,2,C,3);
//   for(i=0; i<9; i++)
//     printf("%lf ", C[i]);
//   printf("\n");
// }





knnresult knn_v0(double * X, double * Y, int n, int m, int d, int k)
{

    bool __show_info = true;

    // Start Timer
    time_t t;
    srand((unsigned) time(&t));
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);



    double * XT;    // X Transpose
    mat_transpose(&X, &XT, n, d);

    double * YT;    // Y Transpose
    mat_transpose(&Y, &YT, m, d);



    /**
     * dgemm:  C = α . A x B + β . C
     * Transpose: CblasNoTrans / CblasTrans
     **/

    int mm, nn, kk, alpha, beta;



    /******************************
     *       D2 = X * Y.'   
     *          (NxM)
     ******************************/

    double * D2 = calloc(m*n, sizeof(double));


    printf("1\n");

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
        n               // N
    );

    printf("dd %d\n", n*d);
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



    /**************************
     *       X .* X
     *        (NxD)
     **************************/
    double* X_X = malloc(n*d*sizeof(double));
    if(X_X==NULL)
    {
        printf("Failed allocating memory [X_X] (knn_v0).\n");
        exit(EXIT_FAILURE);
    }

    
    printf("2\n");
    return;

    for(int i=0; i<n*d; i++)
    {
        X_X[i] = X[i] * X[i];
    }

    

    /**************************
     *       (Y .* Y).'
     *        (DxM)
     **************************/

    double * Y_Y = calloc(m*d, sizeof(double));

    for(int i=0; i<m*d; i++)
    {
        Y_Y[i] = Y[i] * Y[i];
    }

    double * Y_Y_T = calloc(d*m, sizeof(double));

    mat_transpose(&Y_Y, &Y_Y_T, m, d);


    printf("3\n");


    /***************************************
     *           D1   (Nx1)
     * for i=1:n
     *   temp = 0;
     *     for j=1:d
     *       temp = temp + XX(i,j);
     *   D1(i,1) = temp; 
     **************************************/
    double i_tmp;
    double * D1 = calloc(n,sizeof(double));
    for(int i=0; i<n; i++)
    {
        i_tmp = 0.0;
        for(int j=0; j<d; j++)
        {
            i_tmp += (double) X_X[d*i+j];  // tmp += X(i,j)
        }
        D1[i] = i_tmp;
    }



    printf("4\n");

    /***************************************
     *           D3   (1xM)
     *  for i=1:m
     *    temp = 0;
     *    for j=1:d
     *      temp = temp + YY(j,i);
     *    D3(1,i) = temp; 
     **************************************/

    double * D3 = calloc(m, sizeof(double));
    for(int i=0; i<m; i++)
    {
        i_tmp = 0.0;
        for(int j=0; j<d; j++)
        {
            i_tmp += (double) Y_Y_T[j*m+i];  // tmp += Y_Y_T(i,j)
        }
        D3[i] = i_tmp;
    }




    // printf("--- (Y.*Y).' ---\n");

    // for(int i=0; i<m; i++)
    // {
    //     // printf("%d: ( ", i);
    //     for(int j=0; j<d; j++)
    //     {
            
    //         printf("%d ", (int)mat_read_ij(&Y, i, j, n) );

    //     }

    //     printf("; \n");
        
    // }

    // printf("\n");




    // printf("--- D3 ---\n");
    // for(int i=0; i<m; i++){
    //     printf("%d ", (int)D3[i]);
    // }
    // printf("\n\n");


    // printf("--- D1 ---\n");
    // for(int i=0; i<n; i++){
    //     printf("%d ", (int)D1[i]);
    // }
    // printf("\n\n");



    /***************************
     *     D2 <-- D12  (NxM)
     *      (Nx1) + (NxM)
     ***************************/

    for(int i=0; i<n; i++)
    {
        i_tmp = D1[i];
        for(int j=0; j<m; j++)
        {
            D2[i*m +j] += i_tmp;    // D2(i,j) += D1(i,1)
        }
    }


    /***************************
     *     D2 <-- D23  (NxM)
     *      (1xN) + (NxM)
     ***************************/
    for(int j=0; j<m; j++)
    {
        i_tmp = D3[j];
        for(int i=0; i<n; i++)
        {
            D2[i*m +j] += i_tmp;    // D2(i,j) += D3(1,i)
        }
    }



    /*************************
     *     sqrt(D2)   (NxM)
     *************************/
    for(int i=0; i<n*m; i++)
    {
        D2[i] = (double)sqrt(D2[i]);
    }

    free(D1);
    free(D3);
    free(X_X);
    free(Y_Y);
    free(Y_Y_T);

    //ouble * C = D2;

    // printf("--- C ---\n");

    // for(int i=0; i<n; i++)
    // {
    //     for(int j=0; j<m; j++)
    //     {

    //         printf("%f ", (double) mat_read_ij(&C, i, j, m) );

    //     }

    //     printf("; \n");
        
    // }

    // printf("\n");



    /**
     * C is ready
     **/










    knnresult res;


    // Stop Timer
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    float delta_us = (float) ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000)/ (1000000);
    if(__show_info)
        printf(" > V0 took %f s\n", delta_us);


    return res;
}



