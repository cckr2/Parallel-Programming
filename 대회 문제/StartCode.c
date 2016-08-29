#include <stdio.h>
#include<stdlib.h> 
#include <mpi.h>

#define ROOT 0

int main(int argc, char *argv[])
{
	int i, j, k;
	int nRank, nProcs;
	MPI_Init(&argc, &argv); 

	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	MPI_Request req;
	MPI_Status status;

	int lastR, start, end, remain, jobs;
	int *displs = (int*)malloc(sizeof(int)*nProcs);
	int *count = (int*)malloc(sizeof(int)*nProcs);
	int mycount;

	lastR = nProcs - 1;
	jobs = JOB / nProcs; //JOB ¹Ù²Ù»ï
	remain = JOB % nProcs; //JOB ¹Ù²Ù»ï
	start = nRank * jobs;

	if (nRank == lastR) {
		end = JOB;//JOB ¹Ù²Ù»ï
	}
	else {
		end = start + jobs;
	}

	for (i = 0; i < nProcs; i++) {
		if (i == lastR) {
			count[i] = jobs + remain;
			displs[i] = jobs*i + remain;
		}
		else {
			count[i] = jobs;
			displs[i] = jobs*i;
		}
	}
	
	mycount = count[nRank];

	MPI_Finalize();
	return 0;
}
