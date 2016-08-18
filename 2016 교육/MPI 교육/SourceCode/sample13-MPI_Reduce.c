#include <stdio.h>
#include <mpi.h>
int main(int argc, char *argv[])
{
	int i, nrank, start, end, ROOT = 0;
	double a[9], sum, tsum;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &nrank);

	start = nrank * 3;
	end = start + 2;
	for (i=start; i<end+1; i++) {
		a[i] = i + 1;
		if (i == start)
			printf("rank (%d) ", nrank);
		printf("a[%d] = %.2f ", i, a[i]);
	}
	
	sum = 0.0;
	for (i=start; i<end+1; i++) 
		sum = sum + a[i];

	MPI_Reduce(&sum, &tsum, 1, MPI_DOUBLE, MPI_SUM,ROOT, MPI_COMM_WORLD);

	if (nrank == ROOT)
		printf("\nrank(%d):sum= %.2f.\n", nrank, tsum);
	printf("\n");
	
	MPI_Finalize();
	return 0;
}
