#include <stdio.h>
#include<math.h>
#include <mpi.h>
#define NodeNum 4
#define ROOT 0

#define L 1.0
#define LAMBDA 0.0003
#define PI 3.141592

int main(int argc, char* argv[]){

	int nRank, nProcs;
	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank); 
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs); 
	MPI_Request req;
	MPI_Status status;

	int JOB = L/LAMBDA;
	int lastR,start,end,remain, jobs;

	lastR = NodeNum - 1;
	jobs = JOB / NodeNum; 
	remain = JOB % NodeNum; 
	start = nRank * jobs;
	if (nRank == lastR) {
		end = start + jobs + remain;
	}
	else {
		end = start +  jobs;
	}
	
	long int count = 0, num_modes;
	double x, y, z, r, delta, sol;

	delta = LAMBDA;
	int i;
	if(nRank==lastR){
		for (i= start; i<=end; i++){
			x = i*delta;
			for (y = 0.0; y <= L; y = y + delta){
				for (z = 0.0; z <= L; z = z + delta){
					r = sqrt(x*x + y*y + z*z);
					if (r <= L)
						count++;
				}
			}
		}	
	}
	else{
		for (i = start; i< end; i++){
			x = i*delta;
			for (y = 0.0; y <= L; y = y + delta){
				for (z = 0.0; z <= L; z = z + delta){
					r = sqrt(x*x + y*y + z*z);
					if (r <= L)
						count++;
				}
			}
		}	
	}
	printf("rank : %d,count : %ld\n",nRank,count);
	long int sum = 0.0;
	MPI_Reduce(&count, &sum, 1, MPI_LONG, MPI_SUM, ROOT, MPI_COMM_WORLD);	
	
	if(nRank==ROOT){
		num_modes = 16 * sum;
		sol = 8.0*PI*L *L*L / 3.0 / LAMBDA / LAMBDA / LAMBDA;
		printf("Number of modes = %ld (from mode counting)\n", num_modes);
		printf("Number of modes = %f (from solution)\n", sol);
	}
	MPI_Finalize(); // MPI End
	return 0;
}
