 #!/bin/bash

#export OPENBLAS_NUM_THREADS=8
#export GOTO_NUM_THREADS=8
#export OMP_NUM_THREADS=8
export CILK_NWORKERS=8

# sudo apt update -y && sudo apt install -y glibc-source

export KNN_PRINT=1
export TIMER_PRINT=1
#export MPI_NWORKERS=1

clear && make clean && make all && clear && mpirun -np 3 main.o $1
# clear && make clean && make all && clear && ./main.o $1
