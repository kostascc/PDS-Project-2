#include <stdio.h>
#include <stdlib.h>
#include "mmarket.h"
#include <stdbool.h>
#include <string.h>
#include "mmarket.h"



int main(int argc, char** argv)
{

    /** 
    * Available Threads.
    * Defaults to 8.
    **/
    int __threads = 8;


    char _x_filename[1024];

    bool _x_value_included = true;

    char _y_filename[1024];

    bool _y_value_included = true;



    /** 
    * For each argument given,
    * find the respective configuration
    * parameter and change it
    **/
    for(int i=0; i<argc; i++){


        int _tmp_int = NULL;   // Temporary placeholder
        char _tmp_str[1024];



        /**
         * " -t_ "
         */
        if(sscanf(argv[i], "-t%d", &_tmp_int))
        {
            if(_tmp_int>0 && _tmp_int<64)
            __threads = _tmp_int;
        }


        /**
         * " -x< Path > "
         */
        if(sscanf(argv[i], "-x%s", _tmp_str))
        {
            strcpy(_x_filename, _tmp_str);
        }

        /**
         * " -x < Path > "
         */
        if(strcmp(argv[i],"-x")==0)
        {
            if(i<argc)
                strcpy(_x_filename, argv[i+1]);
        }



        /**
         * " -y< Path > "
         */
        if(sscanf(argv[i], "-y%s", _tmp_str))
        {
            strcpy(_y_filename, _tmp_str);
        }

        /**
         * " -y < Path > "
         */
        if(strcmp(argv[i],"-y")==0)
        {
            if(i<argc)
                strcpy(_y_filename, argv[i+1]);
        }

    }



    /**
    * Configure threads on system runtime level
    **/

    // Cilk:
    char _tmp_str[4];
    sprintf(_tmp_str, "%d", __threads);
    __cilkrts_set_param("nworkers", _tmp_str);

    // OpenMP:   
    omp_set_num_threads(__threads);



    


    /**
     * Corpus Points
     * n-by-d
     **/
    double * X; 

    /**
     * Query Points
     * m-by-d
     **/
    double * Y;

    /**
     * Corpus points count
     **/
    int n;

    /**
     * Query points count
     **/
    int m;

    /**
     * Dimension count
     * 
     * ! It has to be initialized as 0,
     * otherwise mmarket_import might keep an 
     * arbitrary number of dimensions.
     **/
    int d = 0;
    
    /**
     * neighbours count
     **/
    int k;


    /******************
     ** Create Matrix **
    ******************/

    mmarket_import(_x_filename, &X, &n, &d, _x_value_included, true); // Import MM

    printf("n: %d, d: %d\n", n, d);

    mmarket_import(_y_filename, &Y, &m, &d, _y_value_included, true); // Import MM

    printf("m: %d, d: %d\n", m, d);


    for(int i=0; i<d; i++)
    {
        for(int j=0; j<n; j++)
        {

            printf("%d ", (int)mat_read_ij(&X, i, j, n) );

        }

        printf("\n");
        
    }


printf("\n\n");

    for(int i=0; i<d; i++)
    {
        for(int j=0; j<m; j++)
        {

            printf("%d ", (int)mat_read_ij(&Y, i, j, m) );

        }

        printf("\n");
        
    }


    /*****************
     ** Run versions **
    *****************/

    for(int i=0; i<argc; i++){

        // if(strcmp(argv[i],"-v0")==0){
        //     knn_v0( X , Y , n , m , d , k );
        // }

        

        // if(strcmp(argv[i],"-coorow")==0){
        //     _coo_row(mat);
        // }

        // if(strcmp(argv[i],"-coocol")==0){
        //     _coo_col(mat);
        // }

        // if(strcmp(argv[i],"-coo")==0){
        //     _coo(mat);
        // }
    }



    return 0;

}