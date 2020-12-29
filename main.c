#include <stdio.h>
#include <stdlib.h>
#include "mmarket.h"
#include <stdbool.h>
#include <string.h>
#include "mmarket.h"
#include "knn_v0.h"
#include "knn_v1.h"
#include "knn_v2.h"
#include <mpi.h>
#include "auxlib.h"
#include "mpi_wrapper.h"


// Show Created or imported
// matrices in stdout
#define V0_USE_X_AS_Y true


int main(int argc, char** argv)
{

    /********************
     ** Startup Script **
     ********************/

    _runtime r = startup(argc, argv);



    /******************
     ** Run versions **
     ******************/

    for(int i=0; i<argc; i++){

        if(strcmp(argv[i],"-v0")==0){
            if(V0_USE_X_AS_Y)
            {
                kNN( r.X , r.X , r.n , r.n , r.d , r.k );
            }else{
                kNN( r.X , r.Y , r.n , r.m , r.d , r.k );
            }
        }

        if(strcmp(argv[i],"-v1")==0){
            distrAllkNN( r.X , r.n , r.d , r.k );
        }

        if(strcmp(argv[i],"-v2")==0){
            distrAllkNNVPT( r.X , r.n , r.d , r.k );
        }

    }


    /*************
     ** Cleanup **
     *************/
    
    free(r.X);
    free(r.Y);
    
    if( _MODE_V1_RUNNING || _MODE_V2_RUNNING )
    {
        // Stop MPI Node
        mpi_finish_local();

        // It shouldn't get to this
        return EXIT_FAILURE;
    }

    return 0;
    

}