FC=gfortran
MFC=mpif90
CC=gcc
MCC=mpicc
FFLAGS=-O3
CFLAGS=-O3 

.PHONY: all

all: cserial cparallel fserial fparallel

cparallel: parallel.c
	$(MCC) $(CFLAGS) parallel.c -std=c99 -o c_parallel.exe
cserial: sequential.c
	$(CC) $(CFLAGS) sequential.c -lm -std=c99 -o c_serial.exe
fparallel: parallel.f90
	$(MFC) $(FFLAGS) parallel.f90 -o f_parallel.exe
fserial: sequential.f90
	$(FC) $(FFLAGS) sequential.f90 -o f_serial.exe
clean:
	rm -rf *.o *.mod ?_*l.exe
