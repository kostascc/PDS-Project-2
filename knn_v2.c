#include "knn_v2.h"


knnresult distrAllkNNVPT(double * X, int n_all, int d, int k)
{

    // Start Timer
    time_t t;
    srand((unsigned) time(&t));
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);



    // Stop Timer
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    float delta_us = (float) ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000)/ (1000000);
    

    knnresult knn;
    return knn;

}
