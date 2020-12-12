#include "mat.h"


void mat_set_ij(double** X, double val, int i, int j, int n)
{

    *( *(X) + j + n * i ) = (double)val;

}




double mat_read_ij(double** X, int i, int j, int n)
{

    return (double) *( *(X) + j + n * i ) ;

}


