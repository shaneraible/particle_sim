#
CC = g++
MPCC = mpic++
#Note: this is the flag for gnu compilers. Change this to -openmp for Intel compilers
OPENMP = -fopenmp
MPI = 
CFLAGS = 
LIBS = -lm
NVCC = nvcc
NVCCFLAGS = -O3 -arch=compute_35 -code=sm_35 -gencode arch=compute_35,code=[compute_35,sm_35] -gencode arch=compute_61,code=[compute_61,sm_61]


TARGETS = serial openmp autograder

all:	$(TARGETS)

gpu: gpu.o common.o bin_grid.o
	$(NVCC) -o $@ $(NVCCLIBS) gpu.o common.o bin_grid.o 

serial: serial.o common.o bin_grid.o
	$(CC) -o $@ serial.o common.o bin_grid.o -lm
autograder: autograder.o common.o
	$(CC) -o $@ autograder.o common.o -lm
openmp: openmp.o common.o
	$(CC) -o $@ $(OPENMP) openmp.o common.o -lm
mpi: mpi.o common.o
	$(MPCC) -o $@ $(MPILIBS) mpi.o common.o -lm

gpu.o: gpu.cu common.h
	$(CC) -c $(NVCCFLAGS) gpu.cu

autograder.o: autograder.cpp common.h
	$(CC) -c $(CFLAGS) autograder.cpp
openmp.o: openmp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp
mpi.o: mpi.cpp common.h
	$(MPCC) -c $(MPI) $(CFLAGS) mpi.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp
bin_grid.o: bin_grid.cpp bin_grid.h
	$(CC) -c $(CFLAGS) bin_grid.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
