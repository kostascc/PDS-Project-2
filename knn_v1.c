/**
 * PDS Project 2
 * 
 * Copyright 2021 Ⓒ K. Chatzis
 * kachatzis <at> ece.auth.gr
 */

#include "knn_v1.h"



knnresult distrAllkNN(double * X, int n_all, int d, int k)
{

    // Start Timer
    time_t t;
    srand((unsigned) time(&t));
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);


    /**
     * Set V0 in V1 mode. This changes
     * printing and memory management 
     * in V0.
     **/
    _MODE_V1_RUNNING = true;

    // Initialize MPI
    int node_id, cluster_size;
    mpi_initialize(&node_id, &cluster_size);


    if(cluster_size<2 && node_id==0)
    {
        printf("There is no cluster!\n");
        mpi_abort();
    }

    if(k > n_all/cluster_size && node_id==0)
    {
        printf("Illegal K!\n");
        mpi_abort();
    }


    // Who I receive from 
    int node_receive = node_id-1;
    if(node_receive<0)
    {
        node_receive = cluster_size-1;
    }

    // Who I send to
    int node_send = node_id+1;
    if(node_send >= cluster_size)
    {
        node_send = 0;
    }


    // Allocate
    int alloc_size = d*((n_all/cluster_size)+cluster_size);

    double* Y = (double*) malloc(alloc_size*sizeof(double));
    if(Y==NULL) mpi_abort_msg("Malloc Failed (Y)");

    double* Z = (double*) malloc(alloc_size*sizeof(double));
    if(Z==NULL) mpi_abort_msg("Malloc Failed (Z)");

    int n = n_per_node(node_id, cluster_size, n_all);

    // Create an array with results
    // that fits the whole cluster's
    // Query points.
    knnresult* res = (knnresult*)malloc(cluster_size*sizeof(knnresult));

    // Request array for receive and send
    MPI_Request mpi_request[2];


    /**
     * TODO: Give option for local corpus
     * loading instead of MPI Based.
     **/

    double *XX = NULL;
    
    // if(DEBUG_CHK_BATCH)
    // {
    //     XX = malloc(n_all*d*sizeof(double));
    //     if(XX==NULL) mpi_abort_msg("Malloc Failed (XX)");
    //     memcpy(&XX[0], &X[0], n_all*d*sizeof(double));
    // }


    // Nodes other than the main, should
    // receive the X matrix first.
    // Node 0 already has it.
    // if(node_id>0)
    // {
    //     mpi_receive_data_b(MPI_MODE_DATA_DISTRIBUTION, X, n_all*d, node_receive, mpi_request);
    // }

    // All nodes, except the last one, 
    // should send the matrix to the
    // next one.
    // if(node_id != cluster_size-1)
    // {
    //     mpi_send_data_b(MPI_MODE_DATA_DISTRIBUTION, X, n_all*d, node_send, mpi_request);
    // }

    // Hold on to local data only
    memcpy(&X[0], &X[ d*(n_all/cluster_size)*node_id ], d*n_per_node(node_id, cluster_size, n_all) *sizeof(double));
    X = realloc(X, d*n_per_node(node_id, cluster_size, n_all)*sizeof(double));
    if(X==NULL) mpi_abort_msg("Realloc Failed (X)");

    // Copy X to Y
    memcpy(&Y[0], &X[0], d*n_per_node(node_id, cluster_size, n_all) *sizeof(double));


    // Which part of Query
    // Y am I receiving now
    int working_batch = node_id;

    // Current query points on Y
    int m;

    // For each batch of Y query points
    for(int i_rec=0; i_rec<cluster_size; i_rec++)
    {
        
        // if(DEBUG_A==1)
        //     printf("\n\n\nNode %d working on batch %d\n", node_id, working_batch);

        // Wait for the last message
        // to finish receiving.
        // int batch_offset = working_batch*(n_all/cluster_size);

        if(working_batch!=node_id)
        {
            mpi_send_data_wait(mpi_request);
            mpi_receive_data_wait(mpi_request);
            memcpy(&Y[0], &Z[0], m*d*sizeof(double));
        }

        m = n_per_node(working_batch, cluster_size, n_all);


        // Send Current Working Batch
        mpi_send_data_nb(MPI_MODE_CORPUS_DISTRIBUTION, Y, m*d, node_send, mpi_request);


        // wte: What To Expect
        // What length to expect, depending
        // on the next batch id
        int wte = working_batch-1;
        if(wte<0)
            wte = cluster_size -1;
        wte = n_per_node(wte, cluster_size, n_all);

        // Receive Next batch
        // mpi_receive_data_nb(MPI_MODE_CORPUS_DISTRIBUTION, Z, wte, node_receive, mpi_request);
        mpi_receive_data_nb(MPI_MODE_CORPUS_DISTRIBUTION, Z, wte*d, node_receive, mpi_request);

        res[working_batch] = kNN(X, Y, n, m, d, k);
        // res[working_batch].k = k;    // Override k
        // res[working_batch].m = m;    // Override m

        // printf("\n^^ (%d) Batch: %d, m: %d, offset: %d ^^\n", node_id, working_batch, m,  node_id*n_all/cluster_size);
        // Offset index results based on batch id
        // TODO: Make it parallel
        cilk_spawn knnresult_offset_nidx(&(res[working_batch]), node_id*(n_all/cluster_size));


        // update working batch id
        working_batch--;
        if(working_batch<0)
            working_batch = cluster_size-1;

        // update batch's query point count
        m = n_all/cluster_size;
        m += (working_batch==cluster_size-1)? n_all%cluster_size: 0;

        
    }//end: for i=0:cluster_size-1


    knnresult* knn = (knnresult*)malloc(sizeof(knnresult));
    knn->m = n_all;
    knn->k = k;

    cilk_sync;  // knnresult_offset_nidx(...)

    /**
     * Collects KNN results from other nodes, compares and
     * selects the best neighbors, then distributes them
     * to the next node.
     * For master Nodes, this only collects & selects.
     */
    collect_n_propagate_knn(res, knn, node_id, node_receive, node_send, cluster_size);


    // Clean Up
    for(int i=0; i<cluster_size; i++)
    {
        free(res[i].ndist);
        free(res[i].nidx);
    }
    free(res);
    
    if(X!=NULL)
    {
        free(X);
        X = NULL;
    }

    if(Y!=NULL)
    {
        free(Y);
        Y = NULL;
    }

    if(node_id>0)
    {
        mpi_finalize(); // bye
        return *knn;    // ~~~
    }


    /**
     * Result
     **/
    if(_KNN_PRINT)
    {
        print_res(knn); 
    }

    // Stop Timer
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    float delta_us = (float) ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000)/ (1000000);
    if(_TIMER_PRINT)
    {
        printf(" > V1 took %f s [N:%d, D:%d, K:%d]\n", delta_us, n_all, d, k);
    }

    mpi_finalize();
    return *knn;

}
