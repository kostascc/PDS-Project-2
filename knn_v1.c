#include "knn_v1.h"


#define DEBUG_A 1
#define DEBUG_MPI 1


#define MPI_MODE_INITIALIZING 0
#define MPI_MODE_DATA_DISTRIBUTION 1
#define MPI_MODE_KNN_COLLECTION_DIST 2
#define MPI_MODE_KNN_COLLECTION_INDX 3
#define MPI_MODE_KNN_COLLECTION_M 7
#define MPI_MODE_CORPUS_DISTRIBUTION 6
#define MPI_MODE_CALCULATING 4
#define MPI_MODE_COLLECTING 5


int node_id_tmp=-1;



int _v1_n_per_node(int node_id, int cluster_size, int n)
{

    int ret = (n/cluster_size);
    if(node_id == cluster_size-1)
    {
        return ret + n%cluster_size;
    }
    return ret;
}



void _v1_compare_knnresult(knnresult* res, knnresult* res_)
{

    (res)->k = (res_)->k;
    (res)->m = (res_)->m;


    // TODO: result comparison
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
        printf("%d Sending to %d [%d]\n", node_id_tmp, partner, m);

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
        printf("%d Receiving from %d [%d]\n", node_id_tmp, partner, m);

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


    if(cluster_size<2)
    {
        printf("There is no cluster!\n");
        exit(EXIT_FAILURE);
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


    node_id_tmp = node_id;


    // Allocate D
    int alloc_size = (n_all/cluster_size)+cluster_size;
    double* Y = (double*) malloc(alloc_size*sizeof(double));
    double* Z = (double*) malloc(alloc_size*sizeof(double));
    int n = _v1_n_per_node(node_id, cluster_size, n_all);

    // Create an array with results
    // that fits the whole cluster's
    // Query points.
    knnresult res[cluster_size];


    /**
     * Construct kNN result struct
     * for OMPI
     **/
    // const int nitems = 4;
    // int blocklengths[4] = {1, 1, 1, 1};
    // MPI_Datatype types[4] = {MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT};
    // MPI_Datatype MPI_KNNRESULT;
    // MPI_Aint     offsets[4];
    // offsets[0] = offsetof(knnresult, nidx);
    // offsets[1] = offsetof(knnresult, ndist);
    // offsets[2] = offsetof(knnresult, m);
    // offsets[3] = offsetof(knnresult, k);

    // Request array for receive and send
    MPI_Request mpi_request[2];



    // Nodes other than the main, should
    // receive the X matrix first.
    // Node 0 already has it.
    if(node_id>0)
        _v1_receive_data_b(MPI_MODE_DATA_DISTRIBUTION, X, n_all, node_receive, mpi_request);


    // All nodes, except the last one, 
    // should send the matrix to the
    // next one.
    if(node_id != cluster_size-1)
        _v1_send_data_b(MPI_MODE_DATA_DISTRIBUTION, X, n_all, node_send, mpi_request);
    

    // Hold on to local data only
    memcpy(&X[0], &X[ (n_all/cluster_size)*node_id ], _v1_n_per_node(node_id, cluster_size, n_all) *sizeof(double));

    // Xopy X to Y
    memcpy(&Y[0], &X[0], _v1_n_per_node(node_id, cluster_size, n_all) *sizeof(double));


    if(node_id==0)
    {
        for(int i=0; i<cluster_size; i++)
        {
            printf("node %d n: %d\n", i, _v1_n_per_node(i, cluster_size, n_all));
        }
    }

    // Which part of Query
    // Y am I receiving now
    int working_batch_id = node_id;

    int m;


    // For each batch of Y query points
    for(int i_rec=0; i_rec<cluster_size; i_rec++)
    {
        
        
        if(DEBUG_A==1)
            printf("Node %d working on batch %d\n", node_id, working_batch_id);


        
        
        // Wait for the last message
        // to finish receiving.
        if(working_batch_id!=node_id)
            _v1_send_data_wait(mpi_request);


        m = _v1_n_per_node(working_batch_id, cluster_size, n_all);


        // Send Current Working Batch
        _v1_send_data_nb(MPI_MODE_CORPUS_DISTRIBUTION, Y, m, node_send, mpi_request);
        

        // wte: What To Expect
        // What data to expect, depending
        // on the next batch id
        int wte = working_batch_id-1;
        if(wte<0)
            wte = cluster_size -1;
        wte = _v1_n_per_node(wte, cluster_size, n_all);

        // Receive Next batch
        _v1_receive_data_nb(MPI_MODE_CORPUS_DISTRIBUTION, Z, wte, node_receive, mpi_request);



        res[working_batch_id] = kNN(X, Y, n, m, d, k);


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


        

        // }
        // else
        // {
        //     // Non-Blocking Receive & Send
        //     _v1_receive_query( &Z );
        //     _v1_send_query( Y );

        //     // This is non blocking

        // }

        
    }



    // Offset results based on node id
    // TODO: ^


    // Free Up Memory
    // free(X);
    // free(Y);
    // free(Z);

    _v1_send_data_wait(mpi_request);
    _v1_receive_data_wait(mpi_request);
    printf("--------------\n");

    knnresult res_buff[cluster_size];



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
            
            int m_buff[1];
            // Receive M
            _v1_receive_data_b_t(
                MPI_MODE_KNN_COLLECTION_M, 
                m_buff, 
                1, 
                node_receive, 
                MPI_INT, 
                mpi_request
            );

            printf("[%d] batch: %d, m: %d\n", node_id, i, m_buff[0]);

            if(m_buff[0] < 0)
            {
                printf("M received from node %d on node %d was illegal!\n");
                exit(EXIT_FAILURE);
            }

            // TODO: Maybe execute 2 concurrent receives
                printf("indx\n");
            int * indx_buff = (int*)malloc(m_buff[0]*sizeof(int));
            _v1_receive_data_b_t(
                MPI_MODE_KNN_COLLECTION_INDX, 
                indx_buff,
                m_buff[0],
                node_receive,
                MPI_INT,
                mpi_request
            );

            double * dist_buff = (double*)malloc(m_buff[0]*sizeof(double));
            _v1_receive_data_b_t(
                MPI_MODE_KNN_COLLECTION_DIST, 
                dist_buff,
                m_buff[0],
                node_receive,
                MPI_DOUBLE,
                mpi_request
            );

            cilk_sync;
            
            knnresult res_tmp;
            res_tmp.ndist = dist_buff;
            res_tmp.nidx  = indx_buff;
            res_tmp.m     = m_buff[0];
            res_tmp.k     = k;

            cilk_spawn _v1_compare_knnresult(&(res[i]), &(res_tmp));

        }
            
    }//end: if(node_id>1)
        
    // nodes that send
    if(node_id!=0)
    {
        // TODO: Move sending right after calculation
        // to make it async.

        // Send result to next node
        for(int i=0; i<cluster_size; i++)
        {
            int m_buff[] = {res[i].m};
            printf("send m: %d\n", res[i].m);
            _v1_send_data_b_t(
                MPI_MODE_KNN_COLLECTION_M, 
                m_buff,
                1,
                node_send,
                MPI_INT,
                mpi_request
            );
            _v1_send_data_b_t(
                MPI_MODE_KNN_COLLECTION_INDX, 
                res[i].nidx,
                res[i].m,
                node_send,
                MPI_INT,
                mpi_request
            );
            _v1_send_data_b_t(
                MPI_MODE_KNN_COLLECTION_DIST, 
                res[i].ndist,
                res[i].m,
                node_send,
                MPI_DOUBLE,
                mpi_request
            );
        }
    }
        


    // Nodes other than master have finished
    if(node_id>0)
    {
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

        // _v1_send_data_nb(MPI_MODE_KNN_DISTRIBUTION_DIST, res, k*cluster_size, node_send, mpi_request);

        // for(int i=cluster_size-node_id; i<cluster_size-1; i++)
        // {

        //     // receive frrom previous
        //     _v1_receive_data_nb(MPI_MODE_KNN_DISTRIBUTION_DIST, res_buff, k*cluster_size, node_send, mpi_request);
            
        //     // Wait for "res" buffer to be freed
        //     _v1_send_data_wait(mpi_request);

        //     // Copy receive buffer to "res"
        //     memcpy(&res[0], &res_buff[0], cluster_size*sizeof(knnresult));

        //     // send "res" to next
        //     _v1_send_data_nb(MPI_MODE_KNN_DISTRIBUTION_DIST, res, k*cluster_size, node_send, mpi_request);
        
        // }
   
    

    // // Calculate the result for
    // // each owns' query points.
    // res[node_id] = kNN( X, // Corpus
    //                     X, // Query
    //                     n, // Corpus Size
    //                     n, // Query Size
    //                     d, // Dimensions
    //                     k  // Neighbors
    // );





    MPI_Finalize();

    // Stop Timer
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    float delta_us = (float) ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000)/ (1000000);
    

        
    if(_TIMER_PRINT)
    {
        printf(" > V1 took %f s\n", delta_us);
    }

    return res[0];

}


