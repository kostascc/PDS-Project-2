#include "knn_v1.h"



void abort()
{
    if(DEBUG_MPI)
        printf("\nMPI Abort (%d)...\n", _v1_node_id_local);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

void finish_local()
{
    if(DEBUG_MPI)
        printf("\nNode Successful Finish (%d)...\n", _v1_node_id_local);
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}


void _v1_check_offset(knnresult res, int min, int max)
{
    int mmin = max;
    int mmax = min;
    for(int i=0; i<res.m*res.k; i++)
    {
        if(res.nidx[i]>mmax)
            mmax = res.nidx[i];

        if(res.nidx[i]<mmin)
            mmin = res.nidx[i];

        if(res.nidx[i]<min || res.nidx[i]>max)
        {
            printf("Offset Check Failed! (indx:%d, min:%d, max:%d, i:%d)\n", res.nidx[i], min, max, i);
            abort();
        }
    }
    printf("mmax: %d (%d)\nmmin: %d (%d)\n", mmin, min, mmax, max);
}

void _v1_check_batch(double *Y, double *X, int m, int offset)
{
    // m = n * d
    for(int i=0; i<m; i++)
    {
        if(Y[i] != X[i+offset])
        {
            printf("Batch Check Failed! (found:%.2f (%.2f), m:%d, offset:%d, i:%d)\n", Y[i], X[i+offset], m, offset, i);
            abort();
        }
    }
}


void _v1_offset_nidx(knnresult* res, int offset)
{
    for(int i=0; i<res->m * res->k; i++)
    {
        res->nidx[i] += offset;
    }
}


int _v1_n_per_node(int node_id, int cluster_size, int n)
{

    int ret = (n/cluster_size);
    if(node_id == cluster_size-1)
    {
        return ret + n%cluster_size;
    }
    return ret;
}


/**
 * Compare thhe results of different nodes on the same
 * Query, and find the resulting knn.
 * res  : Resulting data.
 * res_ : Comparison data.
 **/
void _v1_compare_knnresult(knnresult* res, knnresult* res_)
{

    res->m = (res->m > res_->m)? res->m : res_->m ;

    int alloc = 2*res->k;


    for(int m = 0; m<res->m; m++)
    {

        double * distm 	= (double *) malloc(alloc*sizeof(double));
        int * idxm 		= (int *)    malloc(alloc*sizeof(int));

        for(int i=0; i<res->k; i++)
        {
            distm[i] = (double) res->ndist[m*res->k+i];
            idxm[i]  =    (int) res->nidx[m*res->k+i] ;
        }

        for(int i=0; i<res_->k; i++)
        {
            distm[i+res->k] = (double) res_->ndist[m*res->k+i];
            idxm[i+res->k]  =    (int) res_->nidx[m*res->k+i] ;
        }

        aux_mergeSort(&distm, &idxm, alloc);

        for(int i=0; i<res->k; i++)
        {
            res->ndist[m*res->k+i] = (double) distm[i];
            res->nidx[m*res->k+i]  =    (int) idxm[i] ;
        }

    }


    // TODO: We must compare different position of k

    // for(int i=0; i< res->k * res->m ; i++)
    // {
    //     bool change = false;
    //     if((res->ndist)[i] > (res_->ndist)[i])
    //         change = true;
        
    //     if(change)
    //     {
    //         (res->nidx)[i]  = (int)(res_->nidx)[i];
    //         (res->ndist)[i] = (double)(res_->ndist)[i];
    //     }
    // }
}


/**
 * Usage:
 * 
 * for non blocking:
 * _v1_send_data_nb(MPI_MODE_DISTRIBUTING, X, n, node_send, mpi_request);
 * _v1_send_data_wait(mpi_request); // call later
 * 
 * for blocking:
 * _v1_send_data_b(...)
 **/


void _v1_send_data_nb_t(int code, double* Y, int m, int partner, MPI_Datatype type, MPI_Request request[])
{

    if(DEBUG_MPI)
        printf("%d Sending to %d [%d]\n", _v1_node_id_local, partner, m);

    MPI_Isend(
        Y,         // *Buffer
        m,          // Count
        type,     // Type
        partner,   // Destination
        code,    // Tag
        MPI_COMM_WORLD,   // Comm
        &(request[1])   // *Request
    );
    
}

void _v1_receive_data_nb_t(int code, double* Z, int m, int partner, MPI_Datatype type, MPI_Request request[])
{

    if(DEBUG_MPI)
        printf("%d Receiving from %d [%d]\n", _v1_node_id_local, partner, m);

    MPI_Irecv(
        Z,         // *Buffer
        m,          // Count
        type,     // Type
        partner,   // Destination
        code,    // Tag
        MPI_COMM_WORLD,   // Comm
        &(request[0])   // *Request
    );

}

// With Default DataType: Double
void _v1_send_data_nb(int code, double* Y, int m, int partner, MPI_Request request[])
{

    _v1_send_data_nb_t(code, Y, m, partner, MPI_DOUBLE, request);
    
}

// With Default DataType: Double
void _v1_receive_data_nb(int code, double* Z, int m, int partner, MPI_Request request[])
{

    _v1_receive_data_nb_t(code, Z, m, partner, MPI_DOUBLE, request);

}

void _v1_send_data_wait(MPI_Request request[])
{

    MPI_Status status;
    MPI_Wait(&(request[1]), &status);
    
}

void _v1_receive_data_wait(MPI_Request request[])
{

    MPI_Status status;
    MPI_Wait(&(request[0]), &status);

}

void _v1_send_data_b(int code, double* Y, int m, int partner, MPI_Request request[])
{

    _v1_send_data_nb(code, Y, m, partner, request);
     _v1_send_data_wait(request);
    
}

void _v1_receive_data_b(int code, double* Z, int m, int partner, MPI_Request request[])
{

    _v1_receive_data_nb(code, Z, m, partner, request);
    _v1_receive_data_wait(request);

}


void _v1_send_data_b_t(int code, double* Y, int m, int partner, MPI_Datatype type, MPI_Request request[])
{

    _v1_send_data_nb_t(code, Y, m, partner, type, request);
     _v1_send_data_wait(request);
    
}

void _v1_receive_data_b_t(int code, double* Z, int m, int partner, MPI_Datatype type, MPI_Request request[])
{

    _v1_receive_data_nb_t(code, Z, m, partner, type, request);
    _v1_receive_data_wait(request);

}








knnresult distrAllkNN(double * X, int n_all, int d, int k)
{

    // Start Timer
    time_t t;
    srand((unsigned) time(&t));
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);


    // Set V0 in V1 mode
    _MODE_V1_RUNNING = true;


    // Initiate MPI
    MPI_Init(NULL, NULL);

    // Get Cluster Size
    int cluster_size;
    MPI_Comm_size(MPI_COMM_WORLD, &cluster_size);

    // Get Node ID
    int node_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &node_id);

    // Set Global Node ID for debugging
    _v1_node_id_local = node_id;

    if(cluster_size<2)
    {
        printf("There is no cluster!\n");
        abort();
    }

    if(k > n_all/cluster_size)
    {
        printf("Illegal K!\n");
        abort();
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
    double* Z = (double*) malloc(alloc_size*sizeof(double));
    int n = _v1_n_per_node(node_id, cluster_size, n_all);

    // Create an array with results
    // that fits the whole cluster's
    // Query points.
    knnresult res[cluster_size];

    // Request array for receive and send
    MPI_Request mpi_request[2];


    // char padd[cluster_size*2];
    // padd[0] = '\0';
    // for(int i=0; i<2*node_id; i++)
    // {
    //     padd[i] = ' ';
    //     padd[i+1] = '\0';
    // }


    double *XX = malloc(n_all*d*sizeof(double));
    memcpy(&XX[0], &X[0], n_all*d*sizeof(double));

    // Nodes other than the main, should
    // receive the X matrix first.
    // Node 0 already has it.
    if(node_id>0)
    {
        _v1_receive_data_b(MPI_MODE_DATA_DISTRIBUTION, X, n_all*d, node_receive, mpi_request);
    }

    // All nodes, except the last one, 
    // should send the matrix to the
    // next one.
    if(node_id != cluster_size-1)
    {
        _v1_send_data_b(MPI_MODE_DATA_DISTRIBUTION, X, n_all*d, node_send, mpi_request);
    }

    // Hold on to local data only
    memcpy(&X[0], &X[ d*(n_all/cluster_size)*node_id ], d*_v1_n_per_node(node_id, cluster_size, n_all) *sizeof(double));

    // Copy X to Y
    memcpy(&Y[0], &X[0], d*_v1_n_per_node(node_id, cluster_size, n_all) *sizeof(double));


    // Which part of Query
    // Y am I receiving now
    int working_batch_id = node_id;

    int m;

    // For each batch of Y query points
    for(int i_rec=0; i_rec<cluster_size; i_rec++)
    {
        
        
        // if(DEBUG_A==1)
        //     printf("\n\n\nNode %d working on batch %d\n", node_id, working_batch_id);


        // Wait for the last message
        // to finish receiving.
        if(working_batch_id!=node_id)
            _v1_send_data_wait(mpi_request);


        m = _v1_n_per_node(working_batch_id, cluster_size, n_all);


        // Send Current Working Batch
        _v1_send_data_b(MPI_MODE_CORPUS_DISTRIBUTION, Y, m*d, node_send, mpi_request);


        // wte: What To Expect
        // What data to expect, depending
        // on the next batch id
        int wte = working_batch_id-1;
        if(wte<0)
            wte = cluster_size -1;
        wte = _v1_n_per_node(wte, cluster_size, n_all);

        // Receive Next batch
        // _v1_receive_data_nb(MPI_MODE_CORPUS_DISTRIBUTION, Z, wte, node_receive, mpi_request);
        _v1_receive_data_b(MPI_MODE_CORPUS_DISTRIBUTION, Z, wte*d, node_receive, mpi_request);

        res[working_batch_id] = kNN(X, Y, n, m, d, k);
        res[working_batch_id].k = k;
        res[working_batch_id].m = m;
        // printf("\n^^ (%d) Batch: %d, m: %d, offset: %d ^^\n", node_id, working_batch_id, m,  node_id*n_all/cluster_size);
        // Offset index results based on batch id
        // TODO: Make it parallel
        _v1_offset_nidx(&(res[working_batch_id]), node_id*((int)n_all/cluster_size));

        if(DEBUG_A)
        {
            printf(":: (%d) Batch: %d, m: %d, offset: %d ::\n", node_id, working_batch_id, m,  node_id*n_all/cluster_size);

            printf("   { ");
            for(int kk=0; kk<m*k; kk++)
            {
                printf("%d ", res[working_batch_id].nidx[kk]);
            }
            printf("}\n\n\n");

            print_res( res[working_batch_id] );

            if(DEBUG_CHK_BATCH)
                _v1_check_batch(Y, XX, m, d*working_batch_id*n_all/cluster_size);

            if(DEBUG_CHK_OFFSET)
                _v1_check_offset(res[working_batch_id], node_id*n_all/cluster_size, 
                node_id*n_all/cluster_size+_v1_n_per_node(node_id, cluster_size, n_all) -1);

            printf("^^ (%d) Batch: %d, m: %d, offset: %d ^^\n", node_id, working_batch_id, m,  node_id*n_all/cluster_size);
        }

        /**
         * Update current Batch ID.
         * This show what node the incoming
         * data Z comes from. It goes 
         * anti-clockwise by 1 for each loop.
         **/
        working_batch_id-- ;

        // Make sure it does a cycle
        if(working_batch_id<0)
            working_batch_id = cluster_size-1 ;


        m = _v1_n_per_node(working_batch_id, cluster_size, n_all);


        // wait for sending to finish,
        // so we can change the buffer Y
        _v1_send_data_wait(mpi_request);

        // Copy the incoming batch to
        // the current batch.
        // In the next loop, we will send 
        // and work on this one concurrently,
        // while receiving a new one in Z.
        memcpy(Y, Z, m * sizeof(double));
        
    }//end: for i=0:cluster_size-1

    // cilk_sync; // indx Offset


    // Offset results based on node id
    // TODO: ^


    // Free Up Memory
    // free(X);
    // free(Y);
    // free(Z);

    _v1_send_data_wait(mpi_request);
    _v1_receive_data_wait(mpi_request);
    // if(DEBUG_A)
    //     printf("-------(%d)-------\n", node_id);

    knnresult res_buff[cluster_size];
    
    // free(X);
    // free(XX);
    // free(Y);
    // free(Z);

    // MPI_Finalize();
    // return;


    /**
     * Data Collection:
     * 
     * Worker nodes should receive the
     * kNN result of the previous node
     * (except node=1, doesn't receive)
     * and compare/update with the local
     * results. Next, the receiving nodes
     * should send the result to the next
     * node in the ring.
     * 
     * Master node only receives and compares.
     * All worker nodes have finished
     * after they have sent their kNN
     * results.
     * 
     * The result tranfered is a single 
     * kNN result, therefore a loop is
     * used to account for all batches of
     * Query points.
     **/


    // Nodes that receive
    if(node_id!=1)
    {

        // For each batch
        for(int i=0; i<cluster_size; i++)
        {

            // TODO: Remove m receiving,
            // We already know m is n_all/cluster_size.
            // The last one will have locally the same
            // m plus n_all%cluster_size.
            
            // Receive M
            int m_buff[1];
            _v1_receive_data_b_t(
                MPI_MODE_KNN_COLLECTION_M, 
                m_buff, 
                1, 
                node_receive, 
                MPI_INT, 
                mpi_request
            );

            // if(DEBUG_A)
            //     printf("[%d] batch: %d, m: %d\n", node_id, i, m_buff[0]);

            // if(m_buff[0] < 0)
            // {
            //     printf("M received from node %d on node %d was illegal!\n");
            //     exit(EXIT_FAILURE);
            // }

            // TODO: Maybe execute 2 concurrent receives

            int * indx_buff = (int*)malloc(m_buff[0]*k*sizeof(int));
            _v1_receive_data_b_t(
                MPI_MODE_KNN_COLLECTION_INDX, 
                indx_buff,
                m_buff[0]*k,
                node_receive,
                MPI_INT,
                mpi_request
            );

            double * dist_buff = (double*)malloc(m_buff[0]*k*sizeof(double));
            _v1_receive_data_b_t(
                MPI_MODE_KNN_COLLECTION_DIST, 
                dist_buff,
                m_buff[0]*k,
                node_receive,
                MPI_DOUBLE,
                mpi_request
            );

            knnresult res_tmp;
            res_tmp.ndist = dist_buff;
            res_tmp.nidx  = indx_buff;
            res_tmp.m     = m_buff[0];
            res_tmp.k     = k;

            cilk_spawn _v1_compare_knnresult(&(res[i]), &(res_tmp));

        }

        cilk_sync;
            
    }//end: if(node_id>1)

    // nodes that send
    if(node_id!=0)
    {
        // TODO: Move sending right after calculation
        // to make it async.

        // Send results to next node
        for(int i=0; i<cluster_size; i++)
        {
            // m
            int m_buff[] = { res[i].m };
            _v1_send_data_b_t(
                MPI_MODE_KNN_COLLECTION_M, 
                m_buff,
                1,
                node_send,
                MPI_INT,
                mpi_request
            );

            // nidx
            _v1_send_data_b_t(
                MPI_MODE_KNN_COLLECTION_INDX, 
                res[i].nidx,
                res[i].m*k,
                node_send,
                MPI_INT,
                mpi_request
            );

            // ndist
            _v1_send_data_b_t(
                MPI_MODE_KNN_COLLECTION_DIST, 
                res[i].ndist,
                res[i].m*k,
                node_send,
                MPI_DOUBLE,
                mpi_request
            );

        }
    }



    // Nodes other than master have finished
    if(node_id>0)
    {
        _v1_send_data_wait(mpi_request);
        _v1_receive_data_wait(mpi_request);
        finish_local();
    }


    // Only Master node continues here

    // Concetrate results on a single knnresult

    double* knn_ndist = (double*) malloc(n_all*k*sizeof(double));
    int* knn_nidx     =    (int*) malloc(n_all*k*sizeof(int))   ;

    

    for(int i=0; i<cluster_size; i++)
    {

        int lgth = (i==cluster_size-1) ? 
            (n_all/cluster_size)+n_all%cluster_size : 
            n_all/cluster_size ;

        // printf("\ni %d\n", i);
        // printf("--indxx:\n");
        // for(int j=0; j<res[i].m; j++)
        // {
        //     printf("  ");
        //     for(int p=0; p<res[i].k; p++){
        //         printf("%d, ", res[i].nidx[j*res[i].k + p]);
        //     }
        //     printf("\n");
        // }
        // printf("\n");

        // printf("--distx:\n");
        // for(int j=0; j<res[i].m; j++)
        // {
        //     printf("  ");
        //     for(int p=0; p<res[i].k; p++){
        //         printf("%.2f, ", res[i].ndist[j*res[i].k + p]);
        //     }
        //     printf("\n");

        // }
        // printf("\n");
        // printf("--lgth:%d\n", lgth);

        // Copy batch ndist
        cilk_spawn memcpy(
            &knn_ndist[ i*k*(n_all/cluster_size) ],
            &(res[i].ndist[0]),
            res[i].k * res[i].m * sizeof(double)
        );

        // Copy batch nidx
        cilk_spawn memcpy(
            &knn_nidx[ i*k*(n_all/cluster_size) ],
            &(res[i].nidx[0]),
            res[i].k * res[i].m * sizeof(int)
        );

    }

    cilk_sync; // Sync result concetration

    knnresult knn;
    knn.m = n_all;
    knn.k = k;
    knn.ndist = knn_ndist;
    knn.nidx = knn_nidx;


    /**
     * Result
     **/
    if(_KNN_PRINT)
    {
        
        printf("\n--- RES V1 ---\n");
        print_res(knn);
    }


    MPI_Finalize();

    // Stop Timer
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    float delta_us = (float) ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000)/ (1000000);
    

        
    if(_TIMER_PRINT)
    {
        printf(" > V1 took %f s\n", delta_us);
    }

    return knn;

}


