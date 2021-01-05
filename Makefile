CC=gcc
MPICC=mpicc
CILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/clang
CFLAGS=-O3

default: all

auxlib:
	$(MPICC) $(CFLAGS) -c -fpic auxlib.c -fcilkplus
	$(MPICC) $(CFLAGS) -shared -o libauxlib.so auxlib.o -fcilkplus

knn_v0:
	$(MPICC) $(CFLAGS) -c -fpic knn_v0.c -fcilkplus
	$(MPICC) $(CFLAGS) -shared -o libknn_v0.so knn_v0.o -fcilkplus

knn_v1:
	$(MPICC) $(CFLAGS) -c -fpic knn_v1.c -fcilkplus 
	$(MPICC) $(CFLAGS) -shared -o libknn_v1.so knn_v1.o -fcilkplus

knn_v2:
	$(MPICC) $(CFLAGS) -c -fpic knn_v2.c -fcilkplus 
	$(MPICC) $(CFLAGS) -shared -o libknn_v2.so knn_v2.o -fcilkplus

msort:
	$(CILKCC) $(CFLAGS) -c -fpic msort.c
	$(CILKCC) $(CFLAGS) -shared -o libmsort.so msort.o

mmio:
	$(CC) $(CFLAGS) -c -fpic mmio.c 
	$(CC) $(CFLAGS) -shared -o libmmio.so mmio.o

mmarket:
	$(CC) $(CFLAGS) -c -fpic mmarket.c
	$(CC) $(CFLAGS) -shared -o libmmarket.so mmarket.o

mpi_wrapper:
	$(MPICC) $(CFLAGS) -c -fpic mpi_wrapper.c
	$(MPICC) $(CFLAGS) -shared -o libmpi_wrapper.so mpi_wrapper.o

mat:
	$(CC) $(CFLAGS) -c -fpic mat.c -fcilkplus
	$(CC) $(CFLAGS) -shared -o libmat.so mat.o -fcilkplus

main:
	$(MPICC) $(CFLAGS) -o main.o main.c auxlib.o knn_v0.o  \
	knn_v1.o knn_v2.o mmio.o mmarket.o mat.o mpi_wrapper.o \
	-lopenblas -fcilkplus -fopenmp -lpthread -lm

tester: 
	$(CC) $(CFLAGS) -o tester.o tester.c

all: prog tester

prog: auxlib mpi_wrapper mmio mat mmarket knn_v0 knn_v1 knn_v2 main

.PHONY: all test clean

test:
	./main.o $(ARGS)

clean:
	rm -f *.so
	rm -f *.o
	rm -f *.lib
	rm -f *.txt
