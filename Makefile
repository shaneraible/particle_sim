#
CC = g++
MPCC = mpic++
#Note: this is the flag for gnu compilers. Change this to -openmp for Intel compilers
OPENMP = -fopenmp
MPI = 
CFLAGS = 
LIBS = -lm


TARGETS = serial openmp autograder

all:	$(TARGETS)

serial: serial.o common.o bin_grid.o
	$(CC) -o $@ serial.o common.o bin_grid.o -lm
autograder: autograder.o common.o
	$(CC) -o $@ autograder.o common.o -lm
openmp: openmp.o common.o
	$(CC) -o $@ $(OPENMP) openmp.o common.o -lm
mpi: mpi.o common.o
	$(MPCC) -o $@ $(MPILIBS) mpi.o common.o -lm

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
