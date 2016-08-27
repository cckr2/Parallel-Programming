#include <stdio.h>
#include<time.h> 
#include<stdlib.h> 
#include<math.h>
#include <mpi.h>
int main(int argc, char *argv[])
{
	int nRank, nProcs;
	MPI_Init(&argc, &argv); // MPI Start
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank); // Get current processor rank id
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs); // Get number of processors

	int i, j; 
	float fx[PNUM], fy[PNUM], fz[PNUM], x[PNUM], y[PNUM], z[PNUM];  
	float r3, bufx, bufy, bufz; 
	FILE *fp; 
	fp = fopen("./position.txt", "r"); 
	for (i = 0; i < PNUM; i++){ 
		fscanf(fp, "%f %f %f", &bufx, &bufy, &bufz);
		x[i] = bufx;   y[i] = bufy;   z[i] = bufz;
	}  
	fclose(fp);
	for (i = 0; i < PNUM; i++) {
		fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0;
		for (j = 0; j < PNUM; j++){ 
			if (i == j) 
				continue;  
			r3 = 1.0 / ((x[i] - x[j])*(x[i] - x[j]) + (y[i] - y[j])*(y[i] - y[j]) + (z[i] - z[j])*(z[i] - z[j]));
			r3 *= sqrt(r3);    fx[i] += (x[j] - x[i])*r3;   
			fy[i] += (y[j] - y[i])*r3;    fz[i] += (z[j] - z[i])*r3;
		}
	}  

	for (i = 0; i < PNUM; i++) { 
		printf("%f\t%f\t%f\n", fx[i], fy[i], fz[i]);
	} 

	MPI_Finalize(); // MPI End
	return 0;
}
