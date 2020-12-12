#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <cilk/cilk.h>
#include "mat.h"
#include "mmarket.h"
#include <string.h>
#include <omp.h>
#include "v2.c"
#include "v3_clk.h"
#include "v3_omp.h"
#include "v3.h"
#include "v4.h"
#include "v4_clk.h"
#include "v4_omp.h"
#include "v4_ptd.h"



int main(int argc, char *argv[]){


   /** 
    * Available Threads.
    * Defaults to 8.
    **/
   int __threads = 8;


   /**
    * Prinnt C Vector.
    * Defaults to false.
    **/
   bool __show_c = false;


   /**
    * Print Timing and Status info.
    * Defaults to true.
    **/
   bool __show_info = true;


   /** 
    * For each argument given,
    * find the respective configuration
    * parameter and change it
    **/
   for(int i=0; i<argc; i++){


      int _tmp;   // Temporary placeholder
      

      /**
       * " -t_ "
       * Configures thread count.
       * Valid options are any number 
       * above 1 and below 1024.
       */
      if(sscanf(argv[i], "-t%d", &_tmp))
      {
         if(_tmp>0 && _tmp<1024)
            __threads = _tmp;
      }

      
      /**
       * " -c_ "
       * Configures what data is printed during
       * execution. Valid options are:
       *   0 for info without C vector.
       *   1 for info and C vector.
       *   2 for C vector only.
       * Default Configuration is 0.
       */
      if(sscanf(argv[i], "-c%d", &_tmp))
      {
         if(_tmp==0)
         {
            __show_c = false;
            __show_info = true;
         }
         if(_tmp==1)
         {
            __show_c = true;
            __show_info = true;
         }
         if(_tmp==2)
         {
            __show_c = true;
            __show_info = false;
         }
      }
      
      _tmp = NULL;

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
    * Show Welcome message
    **/
   if(__show_info)
   {
     printf("*************************************************\n");
      printf("*  Triangulinator - Made by K. Chatzis           \n");
      printf("*    <kachatzis@ece.auth.gr>                     \n");
      printf("*  ~ Configuration ~                             \n");
      printf("* Threads: %d                                    \n", __threads);
      printf("*************************************************\n");
   }
   



   /**
    * If the inputs were not enough,
    * show usage information.
    **/
   if (argc < 3)
	{
      printf("-----------------------------------------------\n");
		printf("Usage: %s [martix-market-filename]\n\n", argv[0]);
      printf(" -tX     Set Threads\n\n");
      printf(" -cX     Set display type.\n");
      printf("         X=0: Only Status info (Default).\n");
      printf("         X=1: Status info and C vector.\n");
      printf("         X=2: Only C vector.\n\n");
      printf(" -vX     To run the serialized code.\n");
      printf("         X=2,3,4: For V2, V3, V4 respectively.\n\n");
      printf(" -vXclk  To run OpenCilk code.\n");
      printf("         X=3,4: For V3, V4 respectively.\n\n");
      printf(" -vXomp  To run OpenMP code.\n");
      printf("         X=3,4: For V3, V4 respectively.\n\n");
      printf(" -v4ptd  To run the V4 POSSIX Threads code.\n\n");
      printf(" -coorow Print 1-Based COO rows horizontally.\n");
      printf(" -coocol Print 1-Based COO columns horiizonntally.\n");
      printf(" -coo    Print 1-Based COO rows and columns vertically.\n");
      printf("-----------------------------------------------\n");
      printf("\n\n");

		exit(1);    // Exit
	}




   /******************
   ** Create Matrix **
   ******************/

   int* mat = mmarket_import(argv[1], __show_info); // Import MM



   /*****************
   ** Run versions **
   *****************/

   for(int i=0; i<argc; i++){

      if(strcmp(argv[i],"-v2")==0){
         find_triangles_v2(mat);
      }

      if(strcmp(argv[i],"-v3")==0){
         find_triangles_v3(mat, __show_c, __show_info);
      }

      if(strcmp(argv[i],"-v3clk")==0){
         find_triangles_v3_cilk(mat, __show_c, __show_info, __threads);
      }

      if(strcmp(argv[i],"-v3omp")==0){
         find_triangles_v3_omp(mat, __show_c, __show_info, __threads);
      }

      if(strcmp(argv[i],"-v4")==0){
         v4_simple(mat, __show_c, __show_info);
      }

      if(strcmp(argv[i],"-v4clk")==0){
         v4_cilk(mat, __show_c, __show_info, __threads);
      }

      if(strcmp(argv[i],"-v4omp")==0){
         v4_openmp(mat, __show_c, __show_info, __threads);
      }

      if(strcmp(argv[i],"-v4ptd")==0){
         v4_pthread(mat, __show_c, __show_info, __threads);
      }

      if(strcmp(argv[i],"-coorow")==0){
         _coo_row(mat);
      }

      if(strcmp(argv[i],"-coocol")==0){
         _coo_col(mat);
      }

      if(strcmp(argv[i],"-coo")==0){
         _coo(mat);
      }
   }




   /***************
   ** Clean Up ! **
   ***************/

   free(mat);


   // The End!
   return 0;
}



