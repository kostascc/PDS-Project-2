#include "mpi_wrapper.h"


void mpi_initialize(int* node_id, int* cluster_size)
{
    // Initiate MPI
    MPI_Init(NULL, NULL);

    // Get Cluster Size
    MPI_Comm_size(MPI_COMM_WORLD, cluster_size);

    // Get Node ID
    MPI_Comm_rank(MPI_COMM_WORLD, node_id);

}


void mpi_send_data_nb_t(int code, double* Y, int m, int partner, MPI_Datatype type, MPI_Request request[])
{

    if(DEBUG_MPI)
        printf("Sending to %d [%d]\n", partner, m);

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

void mpi_receive_data_nb_t(int code, double* Z, int m, int partner, MPI_Datatype type, MPI_Request request[])
{

    if(DEBUG_MPI)
        printf("Receiving from %d [%d]\n", partner, m);

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
void mpi_send_data_nb(int code, double* Y, int m, int partner, MPI_Request request[])
{

    mpi_send_data_nb_t(code, Y, m, partner, MPI_DOUBLE, request);
    
}

// With Default DataType: Double
void mpi_receive_data_nb(int code, double* Z, int m, int partner, MPI_Request request[])
{

    mpi_receive_data_nb_t(code, Z, m, partner, MPI_DOUBLE, request);

}

void mpi_send_data_wait(MPI_Request request[])
{

    MPI_Status status;
    MPI_Wait(&(request[1]), &status);
    
}

void mpi_receive_data_wait(MPI_Request request[])
{

    MPI_Status status;
    MPI_Wait(&(request[0]), &status);

}

void mpi_send_data_b(int code, double* Y, int m, int partner, MPI_Request request[])
{

    mpi_send_data_nb(code, Y, m, partner, request);
     mpi_send_data_wait(request);
    
}

void mpi_receive_data_b(int code, double* Z, int m, int partner, MPI_Request request[])
{

    mpi_receive_data_nb(code, Z, m, partner, request);
    mpi_receive_data_wait(request);

}


void mpi_send_data_b_t(int code, double* Y, int m, int partner, MPI_Datatype type, MPI_Request request[])
{

    mpi_send_data_nb_t(code, Y, m, partner, type, request);
     mpi_send_data_wait(request);
    
}

void mpi_receive_data_b_t(int code, double* Z, int m, int partner, MPI_Datatype type, MPI_Request request[])
{

    mpi_receive_data_nb_t(code, Z, m, partner, type, request);
    mpi_receive_data_wait(request);

}

void mpi_abort()
{
    if(DEBUG_MPI)
        printf("\nMPI Abort...\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

void mpi_abort_msg(char* msg)
{
    printf("\nAbort: %s \n", msg);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

void mpi_finish_local()
{
    if(DEBUG_MPI)
        printf("\nNode Successfully Finished...\n");
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}

void mpi_finalize()
{
    if(DEBUG_MPI)
        printf("\nNode Successfully Finished...\n");
    MPI_Finalize();
}