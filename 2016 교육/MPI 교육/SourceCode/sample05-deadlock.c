#include <stdio.h>
#include <mpi.h>
/* if BUF_SIZE > 4KB, deadlock */
#define BUF_SIZE (10)
int main(int argc, char *argv[])
{
	int nprocs, nrank, i, ROOT = 0;
	MPI_Status status;
	double a[BUF_SIZE], b[BUF_SIZE];
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &nrank);
	
	for (i=0; i<BUF_SIZE; i++)
		a[i] = i;
	
	if (nrank == ROOT) {
		for (i=0; i<10; i++)
			printf("before> a[%d] = %.0f, b[%d] = %.0f\n", i, a[i], i, b[i]);
	}

	if (nrank == 0) {
		MPI_Send(a, BUF_SIZE, MPI_DOUBLE, 1, 17, MPI_COMM_WORLD);
		MPI_Recv(b, BUF_SIZE, MPI_DOUBLE, 1, 19, MPI_COMM_WORLD,&status);
	}
	else if (nrank == 1) {
		MPI_Send(a, BUF_SIZE, MPI_DOUBLE, 0, 19, MPI_COMM_WORLD);
		MPI_Recv(b, BUF_SIZE, MPI_DOUBLE, 0, 17, MPI_COMM_WORLD,&status);
	}
	printf("\n\n");
	if (nrank == ROOT) {
		for (i=0; i<10; i++)
			printf("after > a[%d] = %.0f, b[%d] = %.0f\n", i, a[i], i,b[i]);
	}
		
	MPI_Finalize();
	printf("\n");
	return 0;
}
