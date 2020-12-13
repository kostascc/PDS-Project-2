#ifndef mmarket_h__
#define mmarket_h__

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <cilk/cilk.h>
#include "mat.h"
#include "mmio.h"
#include "mmarket.h"
#include <stdbool.h>


/**
 * Import Matrix-Market file
 **/
void mmarket_import(char* filename, double** X, int* n, int* d, bool _x_value_included, bool _x_transpose, bool __show_info);


#endif  // mmarket_h__