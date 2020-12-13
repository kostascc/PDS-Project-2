CC=gcc
MPICC=mpicc
CILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/clang
CFLAGS=-O3

default: all


auxlib:
	$(CC) $(CFLAGS) -c -fpic auxlib.c
	$(CC) $(CFLAGS) -shared -o libauxlib.so auxlib.o

knn_v0:
	$(CC) $(CFLAGS) -c -fpic knn_v0.c -lopenblas -lpthread
	$(CC) $(CFLAGS) -shared -o libknn_v0.so knn_v0.o -lopenblas -lpthread

msort:
	$(CILKCC) $(CFLAGS) -c -fpic msort.c -fcilkplus
	$(CILKCC) $(CFLAGS) -shared -o libmsort.so msort.o

mmio:
	$(CC) $(CFLAGS) -c -fpic mmio.c
	$(CC) $(CFLAGS) -shared -o libmmio.so mmio.o

mmarket:
	$(CC) $(CFLAGS) -c -fpic mmarket.c mat.o
	$(CC) $(CFLAGS) -shared -o libmmarket.so mmarket.o

mat:
	$(CC) $(CFLAGS) -c -fpic mat.c -fcilkplus
	$(CC) $(CFLAGS) -shared -o libmat.so mat.o 


# v4:
# 	$(CC) $(CFLAGS) -c -fpic v4.c
# 	$(CC) $(CFLAGS) -shared -o libv4.so v4.o

# v4_clk:
# 	$(CILKCC) $(CFLAGS) -c -fpic v4_clk.c -fcilkplus -lpthread
# 	$(CILKCC) $(CFLAGS) -shared -o libv4_clk.so v4_clk.o

# v4_omp:
# 	$(CC) $(CFLAGS) -c -fpic v4_omp.c -fopenmp
# 	$(CC) $(CFLAGS) -shared -o libv4_omp.so v4_omp.o

# v4_ptd:
# 	$(CC) $(CFLAGS) -c -fpic v4_ptd.c -lpthread
# 	$(CC) $(CFLAGS) -shared -o libv4_ptd.so v4_ptd.o

# triangles: 
# 	$(CC) $(CFLAGS) -o triangles.o triangles.c mat.o -fcilkplus -fopenmp


# gcc ... -lopenblas -lpthread

main:
	$(CC) $(CFLAGS) -o main.o main.c auxlib.o knn_v0.o mmio.o mmarket.o mat.o -lopenblas -fcilkplus -fopenmp


all: auxlib msort mmio mat mmarket knn_v0 main

.PHONY: all test clean

test:
	./main.o $(ARGS)

clean:
	rm -f *.so
	rm -f *.o
	rm -f *.lib
	rm -f *.txt
