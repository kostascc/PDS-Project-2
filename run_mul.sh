#!/bin/sh

# make clean
# make all

export KNN_PRINT=0
export TIMER_PRINT=1
export DIST_PRINT=0
export CILK_NWORKERS=3

apnd="-rand -d8 -k40 -m1"
fexp="time.exp_"
v0="-v0"
v1="-v1"
v2="-v2"
nmpi=3

# rm -f $fexp

for n in 13000
do
    sleeps=$(expr $n / 400)

    mpirun -np $nmpi ./main.o -n$n $apnd $v2 >> $fexp
    sleep $(echo $sleeps"s")

    mpirun -np $nmpi ./main.o -n$n $apnd $v1 >> $fexp
    sleep $(echo $sleeps"s")

    ./main.o -n$n $apnd $v0 >> $fexp
    sleep $(echo $sleeps"s")
    
    printf "\n" >> $fexp
done

