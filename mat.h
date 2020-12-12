#ifndef mat_h__
#define mat_h__


#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include "msort.h"


/**
 * double** X, double val, int i, int j, int n
 **/
void mat_set_ij(double** X, double val, int i, int j, int n);


/**
 * double** X, int i, int j, int n
 **/
double mat_read_ij(double** X, int i, int j, int n);

#endif  // mat_h__