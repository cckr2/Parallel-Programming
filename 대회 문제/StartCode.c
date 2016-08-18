#include <stdio.h>
#include<stdlib.h> 
#include <mpi.h>
#define NodeNum 8
#define root 0
int main(int argc, char *argv[])
{
	int nRank, nProcs;
	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank); 
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs); 
	MPI_Request req;
	MPI_Status status;

	int lastR, start, end, remain, jobs;
	lastR = NodeNum - 1;
	jobs = JOB / NodeNum; //JOB ¹Ù²Ù»ï
	remain = JOB % NodeNum; //JOB ¹Ù²Ù»ï
	start = nRank * jobs;
	if (nRank = lastR) {
		end = start + jobs + remain;
	}
	else {
		end = start + jobs;
	}
	
	MPI_Finalize();
	return 0;
}
