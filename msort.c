#include "msort.h"
#include <cilk/cilk.h>



/** Make a copy of an array from start to end - 1.
Equivalent to Python's arr[start:end] */
int* slice(int *arr, int start, int end)
{
    int *result = (int *) malloc((end - start) * sizeof(int));
    int i;
    for (i = start; i < end; i++)
    {
        result[i - start] = arr[i];
    }
    return result;
}

/** Merge left and right into result, overwriting everything in result. */
void merge(
    int *I, int *J, /*int *V,*/                     // Result
    int *left_I, int *left_J, /*int *left_V,*/      // Left
    int *right_I, int *right_J, /*int *right_V,*/   // Right
    int leftLen, int rightLen)                      // Length
{


    int i, j;

    // I
    i = 0;
    j = 0;
    while(i < leftLen && j < rightLen)
    {
        if (left_I[i] < right_I[j] ||
            (left_I[i] == right_I[j] && left_J[i] < right_J[j]) 
        )
        {
            I[i + j] = left_I[i];
            J[i + j] = left_J[i];
            // V[i + j] = left_V[i];

            i++;
        }
        else
        {
            I[i + j] = right_I[j];
            J[i + j] = right_J[j];
            // V[i + j] = right_V[j];

            j++;
        }

    }

    for(; i < leftLen; i++)
    {
        I[i + j] = left_I[i];
        J[i + j] = left_J[i];
        // V[i + j] = left_V[i];
    }
    for(; j < rightLen; j++)
    {
        I[i + j] = right_I[j];
        J[i + j] = right_J[j];
        // V[i + j] = right_V[j];
    }



    free(left_I);
    free(right_I);
    free(left_J);
    free(right_J);
    // free(left_V);
    // free(right_V);
}

// in-place merge sort
void mergeSort(int *I, int *J, /*int *V,*/ int len)
{
    if (len <= 1)
    {
        return;
    }
    int *left_I = cilk_spawn slice(I, 0, len / 2 + 1);
    int *right_I = cilk_spawn slice(I, len / 2, len);

    int *left_J = cilk_spawn slice(J, 0, len / 2 + 1);
    int *right_J = cilk_spawn slice(J, len / 2, len);

    // int *left_V = cilk_spawn slice(V, 0, len / 2 + 1);
    // int *right_V = cilk_spawn slice(V, len / 2, len);

    cilk_sync;


    cilk_spawn mergeSort(left_I, left_J, /*left_V,*/ len / 2);
    cilk_spawn mergeSort(right_I, right_J, /*right_V,*/ len - (len / 2));

    cilk_sync;

    merge(I, J, /*V,*/ left_I, left_J, /*left_V,*/ right_I, right_J, /*right_V,*/ len / 2, len - (len / 2));
}



void switcharoo_to_lower_triangle(int *I, int *J, int nz){

    int t;

    for(int i=0; i<nz; i++){
        if(J[i] > I[i])
        {
            t = I[i];
            I[i] = J[i];
            J[i] = t;
        }
    }

}
