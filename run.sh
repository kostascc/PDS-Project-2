 #!/bin/bash

export OPENBLAS_NUM_THREADS=8
export GOTO_NUM_THREADS=8
export OMP_NUM_THREADS=8

clear && make clean && make all && clear && make test ARGS="$1"
