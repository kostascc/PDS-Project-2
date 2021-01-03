#include <stdio.h>   /* printf, stderr, fprintf */
#include <unistd.h>  /* _exit, fork */
#include <stdlib.h>  /* exit */
#include <errno.h>   /* errno */
 


int main(int argc, char **argv)
{

  const char* run_seq = "./main.o \0"; // size: 9
  const char* run_mpi = "mpirun -np 5 main.o \0"; // size: 20

  const char* file_append = " > \0";

  const char* exp_v0 = "v0out.exp_ \0";  // Size: 13
  const char* exp_v1 = "v1out.exp_ \0";  // -
  const char* exp_v2 = "v2out.exp_ \0";  // -

  const char* v0 = " -v0\0";
  const char* v1 = " -v1\0";
  const char* v2 = " -v2\0";

  const char* cmp = "cmp --silent \0";


  putenv("KNN_PRINT=1");
  putenv("TIMER_PRINT=0");



  /*********
   * V0
   *********/
  if (!fork())
  {
    printf("Running: V0\n");

    int offset = 0, size = -1;

    char* cmd = (char*) malloc(300*sizeof(char));

    // Command
    size = strlen(run_seq);
    memcpy(&cmd[offset], run_seq, size*sizeof(char));
    offset += size;
    
    // Arguments
    size = strlen(argv[1]);
    memcpy(&cmd[offset], argv[1], size*sizeof(char));
    offset += size;

    // Version Argument
    size = strlen(v0);
    memcpy(&cmd[offset], v0, size*sizeof(char));
    offset += size;

    // File Append Char
    size = strlen(file_append);
    memcpy(&cmd[offset], file_append, size*sizeof(char));
    offset += size;

    // Export File
    size = strlen(exp_v0);
    memcpy(&cmd[offset], exp_v0, size*sizeof(char));
    offset += size;

    cmd[offset] = '\0';
    
    printf("  > %s\n", cmd);

    int pass = system(cmd);
    free(cmd);
    if(pass==0)
      return;
    
    printf("\nV0 Failed!\n");
    exit(EXIT_FAILURE);
    
    return;
  }else{
    int status;
    wait(&status);  // Wait for fork to join
    if(status!=0)
      exit(EXIT_FAILURE);
  }

  
  
  /*********
   * V1
   *********/
  if (!fork())
  {
    printf("Running: V1\n");

    int offset = 0, size = -1;

    char* cmd = (char*) malloc(300*sizeof(char));

    // Command
    size = strlen(run_mpi);
    memcpy(&cmd[offset], run_mpi, size*sizeof(char));
    offset += size;
    
    // Arguments
    size = strlen(argv[1]);
    memcpy(&cmd[offset], argv[1], size*sizeof(char));
    offset += size;

    // Version Argument
    size = strlen(v1);
    memcpy(&cmd[offset], v1, size*sizeof(char));
    offset += size;

    // File Append Char
    size = strlen(file_append);
    memcpy(&cmd[offset], file_append, size*sizeof(char));
    offset += size;

    // Export File
    size = strlen(exp_v1);
    memcpy(&cmd[offset], exp_v1, size*sizeof(char));
    offset += size;

    cmd[offset] = '\0';
    
    printf("  > %s\n", cmd);

    int pass = system(cmd);
    free(cmd);
    if(pass==0)
      return;
    
    printf("\nV1 Failed!\n");
    exit(EXIT_FAILURE);
    
    return;
  }else{
    int status;
    wait(&status);  // Wait for fork to join
    if(status!=0)
      exit(EXIT_FAILURE);
  }


  /*********
   * V2
   *********/
  if (!fork())
  {
    printf("Running: V2\n");

    int offset = 0, size = -1;

    char* cmd = (char*) malloc(300*sizeof(char));

    // Command
    size = strlen(run_mpi);
    memcpy(&cmd[offset], run_mpi, size*sizeof(char));
    offset += size;
    
    // Arguments
    size = strlen(argv[1]);
    memcpy(&cmd[offset], argv[1], size*sizeof(char));
    offset += size;

    // Version Argument
    size = strlen(v2);
    memcpy(&cmd[offset], v2, size*sizeof(char));
    offset += size;

    // File Append Char
    size = strlen(file_append);
    memcpy(&cmd[offset], file_append, size*sizeof(char));
    offset += size;

    // Export File
    size = strlen(exp_v2);
    memcpy(&cmd[offset], exp_v2, size*sizeof(char));
    offset += size;

    cmd[offset] = '\0';
    
    printf("  > %s\n", cmd);

    int pass = system(cmd);
    free(cmd);
    if(pass==0)
      return;
    
    printf("\nV2 Failed!\n");
    exit(EXIT_FAILURE);
    
    return;
  }else{
    int status;
    wait(&status);  // Wait for fork to join
    if(status!=0)
      exit(EXIT_FAILURE);
  }



  /*********
   * Compare V0 - V1
   *********/
  if (!fork())
  {
    printf("Comparing: V0 to V1\n");

    int offset = 0, size = -1;

    char* cmd = (char*) malloc(300*sizeof(char));

    // Command
    size = strlen(cmp);
    memcpy(&cmd[offset], cmp, size*sizeof(char));
    offset += size;
    
    // Export File
    size = strlen(exp_v0);
    memcpy(&cmd[offset], exp_v0, size*sizeof(char));
    offset += size;

    // Export File
    size = strlen(exp_v1);
    memcpy(&cmd[offset], exp_v1, size*sizeof(char));
    offset += size;

    cmd[offset] = '\0';
    
    printf("  > %s\n", cmd);

    int pass = system(cmd);
    free(cmd);
    if(pass==0)
      return;
    
    printf("\nError on V1\n");
    exit(EXIT_FAILURE);
    
    return;
  }else{
    int status;
    wait(&status);  // Wait for fork to join
    if(status!=0)
      exit(EXIT_FAILURE);
  }


  /*********
   * Compare V0 - V2
   *********/
  if (!fork())
  {
    printf("Comparing: V0 to V2\n");

    int offset = 0, size = -1;

    char* cmd = (char*) malloc(300*sizeof(char));

    // Command
    size = strlen(cmp);
    memcpy(&cmd[offset], cmp, size*sizeof(char));
    offset += size;
    
    // Export File
    size = strlen(exp_v0);
    memcpy(&cmd[offset], exp_v0, size*sizeof(char));
    offset += size;

    // Export File
    size = strlen(exp_v2);
    memcpy(&cmd[offset], exp_v2, size*sizeof(char));
    offset += size;

    cmd[offset] = '\0';
    
    printf("  > %s\n", cmd);

    int pass = system(cmd);
    free(cmd);
    if(pass==0)
      return;
    
    printf("\nError on V2\n");
    exit(EXIT_FAILURE);
    
    return;
  }else{
    int status;
    wait(&status);  // Wait for fork to join
    if(status!=0)
      exit(EXIT_FAILURE);
  }



  printf("\n Tests Successful!\n");


    
//    pid_t  pid;
 
//    pid = fork();
//    if (pid == 0)
//    {
//       /* Θυγατρική διεργασία:
//        * Όταν η fork() επιστρέφει 0, τότε ο εκτελούμενος κώδικας είναι της θυγατρικής διεργασίας
//        *
//        * Μέτρηση ως το δέκα, ανά δευτερόλεπτο
//        */
//       int j;
//       for (j = 0; j < 10; j++)
//       {
//          printf("child: %d\n", j);
//          sleep(2);
//       }
//       _exit(0);  /* Δεν χρησιμοποιείται η exit() αλλά η _exit(0) η οποία δεν καλεί χειριστές εξόδου */
//    }
//    else if (pid > 0)
//    { 
//       /* Μητρική διεργασία:
//        * Όταν η fork() επιστρέφει θετικό ακέραιο, τότε ο εκτελούμενος κώδικας είναι της μητρικής διεργασίας και ο ακέραιος το αναγνωριστικό της μόλις κατασκευασθείσας διεργασίας
//        */
//       int i;
//       for (i = 0; i < 10; i++)
//       {
//          printf("parent: %d\n", i);
//          sleep(1);
//       }
//       exit(0);
//    }
//    else
//    {   
//       /* Σφάλμα
//        * Όταν μία κλήση συστήματος επιστρέφει αρνητικό αριθμό, κάποιο σφάλμα έχει συμβεί
//        * (π.χ. το σύστημα έχει φτάσει στο επιτρεπτό πλήθος διεργασιών).
//        */
//       fprintf(stderr, "can't fork, error %d\n", errno);
//       exit(1);
//    }
}