#include <stdio.h>
#include <mpi.h>
int main(int argc, char *argv[])
{
	int nRank, nProcs;
	char procName[MPI_MAX_PROCESSOR_NAME];
	int nNameLen;

	MPI_Init(&argc, &argv); // MPI Start
	
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank); // Get current processor rank id
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs); // Get number of processors
	MPI_Get_processor_name(procName, &nNameLen);
	
	printf("Hello World. (Process name = %s, nRank = %d, nProcs = %d)\n", procName, nRank, nProcs);
	
	MPI_Finalize(); // MPI End
	return 0;
}
