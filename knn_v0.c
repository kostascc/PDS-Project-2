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

    double * XT;    // X Transpose
    mat_transpose(&X, &XT, n, d);

    double * YT;    // Y Transpose
    mat_transpose(&Y, &YT, m, d);

    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,3,3,2,1,X, 3, Y, 3,2,Y,3);


    printf("\nt\n");

    knnresult res;
    return res;
}



