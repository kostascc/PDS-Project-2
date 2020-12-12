#ifndef mat_h__
#define mat_h__



/** Make a copy of an array from start to end - 1.
Equivalent to Python's arr[start:end] */
extern int* slice(int *arr, int start, int end);

/** Merge left and right into result, overwriting everything in result. */
extern void merge(
    int *I, int *J, /*int *V,*/                     // Result
    int *left_I, int *left_J, /*int *left_V,*/      // Left
    int *right_I, int *right_J, /*int *right_V,*/   // Right
    int leftLen, int rightLen)                      // Length
;

// in-place merge sort
extern void mergeSort(int *I, int *J, /*int *V,*/ int len);

void switcharoo_to_lower_triangle(int *I, int *J, int nz);

#endif