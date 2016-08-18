#include <stdio.h>
#include "mpi.h"
int main (int argc, char* argv[])
{
	/* Initialize the library */
	MPI_Init(&argc, &argv);

	printf("Hello world\n");

	/* Wrap it up. */
	MPI_Finalize();
}
