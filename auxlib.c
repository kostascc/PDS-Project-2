#include "auxlib.h"



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


