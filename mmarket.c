/**
 * PDS Project 2
 * 
 * Copyright 2021 Ⓒ K. Chatzis
 * kachatzis <at> ece.auth.gr
 */

#include "mmarket.h"



void mmarket_import(char* filename, double** X, int* n, int* d, bool _x_value_included, bool _x_transpose, bool __show_info){


    time_t t;
    srand((unsigned) time(&t));
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);


    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;   
    int i;


    /**************
    ** Open File **
    **************/


    if(__show_info)
        printf("[Begin Reading File...]\n");


    if ((f = fopen(filename, "r")) == NULL) {
        printf("Couldn't Open the specified file: \"%s\"\n", filename);
        exit(EXIT_FAILURE);
    }


    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(EXIT_FAILURE);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(EXIT_FAILURE);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(EXIT_FAILURE);
        


    /* reseve memory for matrices */

    // COO Row index
    int* I = (int *) malloc(nz * sizeof(int));
    if(I==NULL) exit(EXIT_FAILURE);

    // COO Column index
    int* J = (int *) malloc(nz * sizeof(int));
    if(J==NULL) exit(EXIT_FAILURE);

    // COO Value
    int* val = (int *) malloc(nz * sizeof(double));
    if(val==NULL) exit(EXIT_FAILURE);


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    int offset = 0;
    int xx, yy;
    int vv;


    if(_x_value_included)
    {
        for (i=0; i<nz; i++)
        {

            fscanf(f, "%d %d %d\n", &xx, &yy, &vv);

            I[i] = xx - 1;
            J[i] = yy - 1;
            val[i] = vv;

        }
    }
    else
    {
        for (i=0; i<nz; i++)
        {

            fscanf(f, "%d %d\n", &xx, &yy);

            I[i] = xx - 1;
            J[i] = yy - 1;
            val[i] = 1;
        }
    }
    



    if (f !=stdin) fclose(f);
    


    /*********************
    ** Create final mat **
    *********************/
    if(__show_info)
        printf("[Begin Creating Matrix...]\n");

    
    


    // Size of memory allocation.
    // This is supposed to help later remove
    // any values that don't belong to the matrix.
    if(!_x_transpose)
    {
        *d = ( *d > 0 ) ?
                *d :            // keep d
                N  ;            // set N

        *n = M;
    }
    else
    {
        *d = ( *d > 0 ) ?
                *d :            // keep d
                M  ;            // set N

        *n = N;
    }


    // memory Allocation Size
    int malloc_size =  (*d) * (*n) ;

    if(malloc_size<1)
    {
        printf("Failed calculating memory allocation (mmarket)!\n");
        exit(EXIT_FAILURE);
    }

    // printf("Malloc: %d\n", malloc_size);


    // Size of X: M * N.
    *X = (double *) calloc( malloc_size , sizeof(double) );
    

    if(X == NULL)
    {
        printf("Failed while allocating memory for X (mmarket).\n"); 
        exit(EXIT_FAILURE); 
    }
    // else if(__show_info)
    // {
        
    //     printf("Matrix of ");

    //     if(malloc_size>=1000000)
    //     {
    //         printf("%.2f MBytes ", ((float)malloc_size/1000000));
    //     }
    //     else if(malloc_size>=1000)
    //     {
    //         printf("%.2f KBytes ", (float)malloc_size/1000);
    //     }
    //     else
    //     {
    //         printf("%d Bytes ", malloc_size);
    //     }       

    //     printf("created ( M: %d, Non-Zeroes: %d).\n", M, nz);
    // }


    int ii, jj, va;

    for(int kk=0; kk<nz; kk++)
    {

        if(!_x_transpose)
        {
            ii = I[kk];
            jj = J[kk];
        }
        else
        {
            ii = J[kk];
            jj = I[kk];
        }

        va = val[kk];
        

        if(ii >= *n)
            continue;

        if(jj >= *d)
            continue;




        *( *X + ii * (*d) + jj ) = (double) va;

    }

    


    // for(int i=0; i<M+2; i++)
    // {
    //     // U
    //     mat[i + 2] = u[i];
    // }
    
    // for(int i=0; i<nz; i++)
    // {
    //     // D
    //     mat[  M + i + 2 + 1 ] = J[i];
    // }



    printf(" > Imported M: %d, N: %d, nz: %d as n: %d, d: %d.\n", M, N, nz, *n, *d);

    /***************
    ** Clean Up ! **
    ***************/

    free(I);
    I = NULL;

    free(J);
    J = NULL;

    free(val);
    val = NULL;


    fflush(stdin);
    fflush(stdout);



    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    float delta_us = (float) ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000)/ (1000000);
    
    if(__show_info)
        printf(" > Import took %f s\n", delta_us);


}


