#include <stdio.h>
#include<time.h> 
#include<stdlib.h> 
#include<math.h>
#include <mpi.h>
#define PNUM 100000
#define NodeNum 8
int main(int argc, char *argv[])
{
	int nRank, nProcs;
	MPI_Init(&argc, &argv); // MPI Start
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank); // Get current processor rank id
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs); // Get number of processors
	MPI_Request req;
	MPI_Status status;

	int i, j; 
	int start, end, remain, jobs;

	jobs = PNUM / NodeNum;
	remain = PNUM % NodeNum;
	start = nRank * jobs;
	if (NodeNum - 1 == nRank) {
		end = start + jobs + remain;
		float temp_fx[jobs + remain], temp_fy[jobs + remain], temp_fz[jobs + remain];
	}
	else {
		end = start + jobs;
		float temp_fx[jobs], temp_fy[jobs], temp_fz[jobs];
	}
	float fx[PNUM], fy[PNUM], fz[PNUM], x[PNUM], y[PNUM], z[PNUM];  
	float r3, bufx, bufy, bufz; 
	

	FILE *fp; 
	fp = fopen("./position.txt", "r"); 

	for (i = 0; i < PNUM; i++){ 
		fscanf(fp, "%f %f %f", &bufx, &bufy, &bufz);
		x[i] = bufx;   y[i] = bufy;   z[i] = bufz;
	}  
	fclose(fp);

	int k = 0;
	for (i = start; i < end; i++) {
		fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0;
		for (j = 0; j < PNUM; j++) {
			if (i == j)
				continue;
			r3 = 1.0 / ((x[i] - x[j])*(x[i] - x[j]) + (y[i] - y[j])*(y[i] - y[j]) + (z[i] - z[j])*(z[i] - z[j]));
			r3 *= sqrt(r3);
			temp_fx[k] += (x[j] - x[i])*r3;
			temp_fy[k] += (y[j] - y[i])*r3;
			temp_fz[k] += (z[j] - z[i])*r3;
		}
		k++;
	}
	if (NodeNum - 1 == nRank) {
		MPI_Isend(&(temp_fx + (jobs)), remain, MPI_FLOAT, 0, 99, MPI_COMM_WORLD, &req);
		MPI_Isend(&(temp_fy + (jobs)), remain, MPI_FLOAT, 0, 99, MPI_COMM_WORLD, &req);
		MPI_Isend(&(temp_fz + (jobs)), remain, MPI_FLOAT, 0, 99, MPI_COMM_WORLD, &req);
	}
	else if (nRank ==0)
	{
		MPI_Recv(&(temp_fx + (PNUM - remain)), remain, MPI_FLOAT, NodeNum - 1, 99, MPI_COMM_WORLD, &status);
		MPI_Recv(&(temp_fy + (PNUM - remain)), remain, MPI_FLOAT, NodeNum - 1, 99, MPI_COMM_WORLD, &status);
		MPI_Recv(&(temp_fz + (PNUM - remain)), remain, MPI_FLOAT, NodeNum - 1, 99, MPI_COMM_WORLD, &status);
	}
	MPI_Gather(&temp_fx, jobs, MPI_FLOAT, fx, jobs, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Gather(&temp_fy, jobs, MPI_FLOAT, fy, jobs, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Gather(&temp_fz, jobs, MPI_FLOAT, fz, jobs, MPI_FLOAT, 0, MPI_COMM_WORLD);


	if (nRank == 0) {
		for (i = 0; i < PNUM; i++) {
			printf("%f\t%f\t%f\n", fx[i], fy[i], fz[i]);
		}
	}

	MPI_Finalize(); // MPI End
	return 0;
}
