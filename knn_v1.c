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


    if(cluster_size<2)
    {
        printf("There is no cluster!\n");
        mpi_abort();
    }

    if(k > n_all/cluster_size)
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


    // char padd[cluster_size*2];
    // padd[0] = '\0';
    // for(int i=0; i<2*node_id; i++)
    // {
    //     padd[i] = ' ';
    //     padd[i+1] = '\0';
    // }


    /**
     * TODO: Give option for local corpus
     * loading instead of MPI Based.
     **/

    double *XX = NULL;
    
    if(DEBUG_CHK_BATCH)
    {
        XX = malloc(n_all*d*sizeof(double));
        if(XX==NULL) mpi_abort_msg("Malloc Failed (XX)");
        memcpy(&XX[0], &X[0], n_all*d*sizeof(double));
    }


    // Nodes other than the main, should
    // receive the X matrix first.
    // Node 0 already has it.
    if(node_id>0)
    {
        mpi_receive_data_b(MPI_MODE_DATA_DISTRIBUTION, X, n_all*d, node_receive, mpi_request);
    }

    // All nodes, except the last one, 
    // should send the matrix to the
    // next one.
    if(node_id != cluster_size-1)
    {
        mpi_send_data_b(MPI_MODE_DATA_DISTRIBUTION, X, n_all*d, node_send, mpi_request);
    }

    // Hold on to local data only
    memcpy(&X[0], &X[ d*(n_all/cluster_size)*node_id ], d*n_per_node(node_id, cluster_size, n_all) *sizeof(double));

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
        /*cilk_spawn*/ knnresult_offset_nidx(&(res[working_batch]), node_id*(n_all/cluster_size));

        // if(DEBUG_A)
        // {
        //     // cilk_sync;  // Joins only if Debugging is required

        //     printf(":: (%d) Batch: %d, m: %d, offset: %d ::\n", node_id, working_batch, m,  node_id*n_all/cluster_size);

        //     printf("   { ");
        //     for(int kk=0; kk<m*k; kk++)
        //     {
        //         printf("%d ", res[working_batch].nidx[kk]);
        //     }
        //     printf("}\n\n\n");

        //     if(DEBUG_RES_I)
        //         print_res( &res[working_batch] );    // TODO: Fix the reference

        // }

        // if(DEBUG_CHK_BATCH)
        // {
        //     printf("  > Batch Offset: %d\n", d*working_batch*(n_all/cluster_size));
        //     knnresult_check_batch(Y, XX, d*m, 
        //         d*working_batch*(n_all/cluster_size)   // Starting Index
        //     );
        // }
            

        // if(DEBUG_CHK_OFFSET)
        // {
        //     // cilk_sync;  // Joins only if Debugging is required
        //     knnresult_check_offset(&res[working_batch], 
        //         node_id*(n_all/cluster_size),     // Minimum
        //         (node_id*(n_all/cluster_size)+n_per_node(node_id, cluster_size, n_all))-1);
        // }
            
        // if(DEBUG_A)
        // {
        //     printf("^^ (%d) Batch: %d, m: %d, offset: %d ^^\n", node_id, working_batch, m,  node_id*n_all/cluster_size);
        // }

        /**
         * Update current Batch ID.
         * This show what node the incoming
         * data Z comes from. It goes 
         * anti-clockwise by 1 for each loop.
         **/
        // working_batch-- ;
        // if(working_batch<0)  // keep it circling
        //     working_batch = cluster_size-1 ;


        // m = n_per_node(working_batch, cluster_size, n_all);


        // wait for sending to finish,
        // so we can change the buffer Y
        // mpi_send_data_wait(mpi_request);

        // Copy the incoming batch to
        // the current batch.
        // In the next loop, we will send 
        // and work on this one concurrently,
        // while receiving a new one in Z.
        // memcpy(&Y[0], &Z[0], d*m * sizeof(double));



        // update working batch id
        working_batch--;
        if(working_batch<0)
            working_batch = cluster_size-1;

        // update batch's query point count
        m = n_all/cluster_size;
        m += (working_batch==cluster_size-1)? n_all%cluster_size: 0;

        
    }//end: for i=0:cluster_size-1

    // cilk_sync; // indx Offset


    // Offset results based on node id
    // TODO: ^


    // Free Up Memory
    // free(X);
    // free(Y);
    // free(Z);

    // cilk_sync;  // Sync Index Offsetting

    // TODO: Not Required?

    // mpi_send_data_wait(mpi_request);
    // mpi_receive_data_wait(mpi_request);
    // // if(DEBUG_A)
    // //     printf("-------(%d)-------\n", node_id);

    // knnresult res_buff[cluster_size];
    
    // // free(X);
    // // free(XX);
    // // free(Y);
    // // free(Z);

    // // MPI_Finalize();
    // // return;


    // /**
    //  * Data Collection:
    //  * 
    //  * Worker nodes should receive the
    //  * kNN result of the previous node
    //  * (except node=1, doesn't receive)
    //  * and compare/update with the local
    //  * results. Next, the receiving nodes
    //  * should send the result to the next
    //  * node in the ring.
    //  * 
    //  * Master node only receives and compares.
    //  * All worker nodes have finished
    //  * after they have sent their kNN
    //  * results.
    //  * 
    //  * The result tranfered is a single 
    //  * kNN result, therefore a loop is
    //  * used to account for all batches of
    //  * Query points.
    //  **/


    // // Nodes that receive
    // if(node_id!=1)
    // {

    //     // For each batch
    //     for(int i=0; i<cluster_size; i++)
    //     {

    //         // TODO: Remove m receiving,
    //         // We already know m is n_all/cluster_size.
    //         // The last one will have locally the same
    //         // m plus n_all%cluster_size.
            
    //         // Receive M
    //         int m_buff[1];
    //         mpi_receive_data_b_t(
    //             MPI_MODE_KNN_COLLECTION_M, 
    //             m_buff, 
    //             1, 
    //             node_receive, 
    //             MPI_INT, 
    //             mpi_request
    //         );

    //         // if(DEBUG_A)
    //         //     printf("[%d] batch: %d, m: %d\n", node_id, i, m_buff[0]);

    //         // if(m_buff[0] < 0)
    //         // {
    //         //     printf("M received from node %d on node %d was illegal!\n");
    //         //     exit(EXIT_FAILURE);
    //         // }

    //         // TODO: Maybe execute 2 concurrent receives

    //         int * indx_buff = (int*)malloc(m_buff[0]*k*sizeof(int));
    //         mpi_receive_data_b_t(
    //             MPI_MODE_KNN_COLLECTION_INDX, 
    //             indx_buff,
    //             m_buff[0]*k,
    //             node_receive,
    //             MPI_INT,
    //             mpi_request
    //         );

    //         double * dist_buff = (double*)malloc(m_buff[0]*k*sizeof(double));
    //         mpi_receive_data_b_t(
    //             MPI_MODE_KNN_COLLECTION_DIST, 
    //             dist_buff,
    //             m_buff[0]*k,
    //             node_receive,
    //             MPI_DOUBLE,
    //             mpi_request
    //         );

    //         knnresult res_tmp;
    //         res_tmp.ndist = dist_buff;
    //         res_tmp.nidx  = indx_buff;
    //         res_tmp.m     = m_buff[0];
    //         res_tmp.k     = k;

    //         cilk_spawn compare_knnresult(&(res[i]), &(res_tmp));

    //     }

    //     cilk_sync;
            
    // }//end: if(node_id>1)

    // // nodes that send
    // if(node_id!=0)
    // {
    //     // TODO: Move sending right after calculation
    //     // to make it async.

    //     // Send results to next node
    //     for(int i=0; i<cluster_size; i++)
    //     {
    //         // m
    //         int m_buff[] = { res[i].m };
    //         mpi_send_data_b_t(
    //             MPI_MODE_KNN_COLLECTION_M, 
    //             m_buff,
    //             1,
    //             node_send,
    //             MPI_INT,
    //             mpi_request
    //         );

    //         // nidx
    //         mpi_send_data_b_t(
    //             MPI_MODE_KNN_COLLECTION_INDX, 
    //             res[i].nidx,
    //             res[i].m*k,
    //             node_send,
    //             MPI_INT,
    //             mpi_request
    //         );

    //         // ndist
    //         mpi_send_data_b_t(
    //             MPI_MODE_KNN_COLLECTION_DIST, 
    //             res[i].ndist,
    //             res[i].m*k,
    //             node_send,
    //             MPI_DOUBLE,
    //             mpi_request
    //         );

    //     }
    // }



    // // Nodes other than master have finished
    // if(node_id>0)
    // {
    //     mpi_send_data_wait(mpi_request);
    //     mpi_receive_data_wait(mpi_request);
    //     mpi_finish_local();
    // }


    // // Only Master node continues here

    // // Concetrate results on a single knnresult

    // double* knn_ndist = (double*) malloc(n_all*k*sizeof(double));
    // int* knn_nidx     =    (int*) malloc(n_all*k*sizeof(int))   ;

    

    // for(int i=0; i<cluster_size; i++)
    // {

    //     int lgth = (i==cluster_size-1) ? 
    //         (n_all/cluster_size)+n_all%cluster_size : 
    //         n_all/cluster_size ;

    //     // printf("\ni %d\n", i);
    //     // printf("--indxx:\n");
    //     // for(int j=0; j<res[i].m; j++)
    //     // {
    //     //     printf("  ");
    //     //     for(int p=0; p<res[i].k; p++){
    //     //         printf("%d, ", res[i].nidx[j*res[i].k + p]);
    //     //     }
    //     //     printf("\n");
    //     // }
    //     // printf("\n");

    //     // printf("--distx:\n");
    //     // for(int j=0; j<res[i].m; j++)
    //     // {
    //     //     printf("  ");
    //     //     for(int p=0; p<res[i].k; p++){
    //     //         printf("%.2f, ", res[i].ndist[j*res[i].k + p]);
    //     //     }
    //     //     printf("\n");

    //     // }
    //     // printf("\n");
    //     // printf("--lgth:%d\n", lgth);

    //     // Copy batch ndist
    //     cilk_spawn memcpy(
    //         &knn_ndist[ i*k*(n_all/cluster_size) ],
    //         &(res[i].ndist[0]),
    //         res[i].k * res[i].m * sizeof(double)
    //     );

    //     // Copy batch nidx
    //     cilk_spawn memcpy(
    //         &knn_nidx[ i*k*(n_all/cluster_size) ],
    //         &(res[i].nidx[0]),
    //         res[i].k * res[i].m * sizeof(int)
    //     );

    // }

    // cilk_sync; // Sync result concetration

    // knnresult knn;
    // knn.m = n_all;
    // knn.k = k;
    // knn.ndist = knn_ndist;
    // knn.nidx = knn_nidx;



    knnresult* knn = (knnresult*)malloc(sizeof(knnresult));
    knn->m = n_all;
    knn->k = k;


    collect_n_propagate_knn(res, knn, node_id, node_receive, node_send, cluster_size);


    /**
     * Worker Nodes should collect knnresults,
     * compare them, and propagate clockwise.
     **/
    if(node_id>0)
    {
        
        // Collect, Compare and Propagate results
        // collect_n_propagate_knn(result, knn, node_id, node_receive, node_send, cluster_size);

        // Clean Up
        for(int i=0; i<cluster_size; i++)
        {
            free(res[i].ndist);
            free(res[i].nidx);
        }
        free(res);

        free(knn);
        
        // End
        // printf("-- END (%d) --", node_id);
        mpi_finish_local(); // bye
    }


    /**
     * Result
     **/
    if(_KNN_PRINT)
    {
        // printf("-- RES V1 --\n");
        print_res(knn);    // TODO: Fix the reference
    }



    // Stop Timer
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    float delta_us = (float) ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000)/ (1000000);
    

        
    if(_TIMER_PRINT)
    {
        printf(" > V1 took %f s\n", delta_us);
    }

    return *knn;

}
