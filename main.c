#include <stdio.h>
#include <stdlib.h>
#include "mmarket.h"
#include <stdbool.h>
#include <string.h>
#include "mmarket.h"
#include "knn_v0.h"
#include "knn_v1.h"
#include <mpi.h>



// Show Created or imported
// matrices in stdout
#define MATRIX_PRINT 0



int main(int argc, char** argv)
{

    // int __threads;

    char _x_filename[1024];

    bool _x_value_included = true;
    bool _x_transpose = false;

    char _y_filename[1024];

    bool _y_value_included = true;
    bool _y_transpose = false;

    int _n = 1000;
    int _m = 100;
    int _d = 10;
    int _k = 5;

    bool _rand = false;



    /** 
    * For each argument given,
    * find the respective configuration
    * parameter and change it
    **/
    for(int i=0; i<argc; i++){


        int _tmp_int = NULL;   // Temporary placeholder
        char _tmp_str[1024];


        /**
         * " -transpx -TraspX -xTransp 
         * -XTransp -xtranspose -xTranspose "
         */
        if( strcmp(argv[i],"-transpx")==0 )
        {
            _x_transpose = true;
            continue;
        }


        /**
         * " -transpy -TraspY -yTransp 
         * -YTransp -ytranspose -yTranspose "
         */
        if( strcmp(argv[i],"-transpy")==0 )
        {
            _y_transpose = true;
            continue;
        }


        /**
         * " -t_ "
         */
        if(sscanf(argv[i], "-t%d", &_tmp_int))
        {
            // if(_tmp_int>0 && _tmp_int<64)
            // {
            //     __threads = _tmp_int;
            //     continue;
            // }
            printf("-t Parameter Deprecated: Please use 'export CILK_NWORKERS=_' instead.\n");
        }



        /**
         * " -k_ "
         */
        if(sscanf(argv[i], "-k%d", &_tmp_int))
        {
            if(_tmp_int>0)
            {
                _k = _tmp_int;
                continue;
            }
        }


        /**
         * " -m_ "
         */
        if(sscanf(argv[i], "-m%d", &_tmp_int))
        {
            if(_tmp_int>0)
            {
                _m = _tmp_int;
                continue;
            }
        }


        /**
         * " -n_ "
         */
        if(sscanf(argv[i], "-n%d", &_tmp_int))
        {
            if(_tmp_int>0)
            {
                _n = _tmp_int;
                continue;
            }
        }


        /**
         * " -d_ "
         */
        if(sscanf(argv[i], "-d%d", &_tmp_int))
        {
            if(_tmp_int>0)
            {
                _d = _tmp_int;
                continue;
            }
        }


        /**
         * " -rand "
         */
        if( strcmp(argv[i],"-rand")==0 )
        {
            _rand = true;
            continue;
        }


        // /**
        //  * " -x< Path > "
        //  */
        // if(sscanf(argv[i], "-x%s", _tmp_str))
        // {
        //     strcpy(_x_filename, _tmp_str);
        //     printf("xfile 1 %s\n", _x_filename);
        //     continue;
        // }


        /**
         * " -x < Path > "
         */
        if(strcmp(argv[i],"-x")==0)
        {
            if(i<argc-1)
            {
                
                strcpy(_x_filename, argv[i+1]);
                continue;
            }
        }


        // /**
        //  * " -y< Path > "
        //  */
        // if(sscanf(argv[i], "-y%s", _tmp_str))
        // {
        //     strcpy(_y_filename, _tmp_str);
        //     printf("yfile 1 %s\n", _y_filename);
        //     continue;
        // }


        /**
         * " -y < Path > "
         */
        if(strcmp(argv[i],"-y")==0)
        {
            if(i<argc-1)
            {
                strcpy(_y_filename, argv[i+1]);
                continue;
            }
        }


    }



    /**
    * Configure threads on system runtime level
    **/

    // Cilk:
    // char _tmp_str[4];
    // sprintf(_tmp_str, "%d", __threads);
    // __cilkrts_set_param("nworkers", _tmp_str);

    // // OpenMP:   
    // omp_set_num_threads(__threads);
    


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
    int n = (_rand) ? _n : 0;

    /**
     * Query points count
     **/
    int m = (_rand) ? _m : 0;

    /**
     * Dimension count
     * 
     * ! It has to be initialized as 0,
     * otherwise mmarket_import might keep an 
     * arbitrary number of dimensions.
     **/
    int d = (_rand) ? _d : 0;
    
    /**
     * neighbours count
     **/
    int k = _k;



    if(!_rand)
    {

        /*******************
         ** Import Matrix **
        ********************/

        mmarket_import(_x_filename, &X, &n, &d, _x_value_included, _x_transpose, true); // Import MM

        printf("n: %d, d: %d\n", n, d);

        mmarket_import(_y_filename, &Y, &m, &d, _y_value_included, _y_transpose, true); // Import MM

        printf("m: %d, d: %d\n", m, d);

    }
    else
    {

        /*******************
        ** Create Matrix **
        ********************/

        // X
        X = (double *) malloc(n*d*sizeof(double));
        if(X==NULL)
        {
            printf("Failed Allocating memory (main).\n");
        }
        for(int i=0; i<n; i++)
        {
            for(int j=0; j<d; j++)
            {
                X[i*d+j] = (double) rand()/RAND_MAX*100.0;
            }
        }


        // Y
        Y = (double *) malloc(n*d*sizeof(double));
        if(Y==NULL)
        {
            printf("Failed Allocating memory (main).\n");
        }
        for(int i=0; i<m; i++)
        {
            for(int j=0; j<d; j++)
            {
                Y[i*d+j] = (double) rand()/RAND_MAX*100.0;
            }
        }

    }


    if(k>n)
    {
        printf("Illegal Parameters (k>n).\n");
        exit(EXIT_FAILURE);
    }

    if(d<1||d>200)
    {
        printf("Illegal Parameters (d<1 or d>200).\n");
        exit(EXIT_FAILURE);
    }

    if(n<1 || m<1)
    {
        printf("Illegal Parameters (n or m <0).\n");
        exit(EXIT_FAILURE);
    }


    if(MATRIX_PRINT)
    {
        printf("\n");

        printf("--- X ---\n");

        for(int i=0; i<n; i++)
        {
            // printf("", i);
            for(int j=0; j<d; j++)
            {

                printf("%d ", (int)mat_read_ij(&X, i, j, d) );

            }

            printf("; \n");
            
        }


        //     mat_transpose(&X, &Y, n, d);

        printf("\n");

        printf("--- Y ---\n");

        for(int i=0; i<m; i++)
        {
            // printf("%d: ( ", i);
            for(int j=0; j<d; j++)
            {

                printf("%d ", (int)mat_read_ij(&Y, i, j, d) );

            }

            printf("; \n");
            
        }

    }


    /*****************
     ** Run versions **
    *****************/

    for(int i=0; i<argc; i++){

        if(strcmp(argv[i],"-v0")==0){
            knn_v0( X , Y , n , m , d , k );
        }

        if(strcmp(argv[i],"-v1")==0){
            knn_v1( X , Y , n , m , d , k );
        }

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


    free(X);
    free(Y);
    

    return 0;

}