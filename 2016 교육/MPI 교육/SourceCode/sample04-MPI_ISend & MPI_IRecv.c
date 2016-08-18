#include <stdio.h>
#include <mpi.h>
int main(int argc, char *argv[])
{
	int nrank, nprocs, tag = 55, ROOT = 0;
	int send = -1, recv = -1;
	MPI_Request req;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &nrank);
	if (nrank == ROOT) {
		printf("Before : nrank(%d) send = %d, recv = %d\n", nrank,send, recv);
		send = 7;
		MPI_Isend(&send, 1, MPI_INTEGER, 1, tag, MPI_COMM_WORLD,&req);
		printf("Other job calculating.\n\n");
		MPI_Wait(&req, &status);
	}
	else {
		MPI_Recv(&recv, 1, MPI_INTEGER, ROOT, tag, MPI_COMM_WORLD,&status);
		printf("After : nrank(%d) send = %d, recv = %d\n", nrank,send, recv);
	}
	MPI_Finalize();
	return 0;
}
