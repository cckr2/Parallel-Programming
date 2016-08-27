FC=gfortran
MFC=mpif90
CC=gcc
MCC=mpicc
FFLAGS=-O3
CFLAGS=-O3 

.PHONY: all

all: cserial cparallel fserial fparallel

cparallel: b.c
	$(MCC) $(CFLAGS) b.c -std=c99 -o c_parallel.exe
cserial: b_serial.c
	$(CC) $(CFLAGS) b_serial.c -lm -std=c99 -o c_serial.exe
fparallel: b.f90
	$(MFC) $(FFLAGS) b.f90 -o f_parallel.exe
fserial: b_serial.f90
	$(FC) $(FFLAGS) b_serial.f90 -o f_serial.exe
clean:
	rm *.o *.mod ?_*l.exe
