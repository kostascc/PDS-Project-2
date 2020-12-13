#include "mat.h"


void mat_set_ij(double** X, double val, int i, int j, int d)
{

    *( *(X) + j + d * i ) = (double)val;

}




double mat_read_ij(double** X, int i, int j, int n)
{

    return (double) *( *(X) + j + n * i ) ;

}


void mat_transpose(double** X, double** T, int M, int N)
{

    *T = malloc(M*N*sizeof(double));

    if(*T==NULL)
    {
        printf("Failed Allocating memory for Transposed Matrix (mat).\n");
        exit(EXIT_FAILURE);
    }


    for(int i=0; i<M; i++)
    {
        for(int j=0; j<N; j++)
        {

            *( *(T) + i + M * j ) = *( *(X) + j + N * i );

        }
    }
}

