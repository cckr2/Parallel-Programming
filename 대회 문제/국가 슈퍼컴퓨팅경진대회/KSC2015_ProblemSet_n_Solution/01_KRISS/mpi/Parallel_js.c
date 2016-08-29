#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <mpi.h>

#define ROOT 0

double genvv(double x){
	return (x*x+pow(x,4)+pow(x,6)+exp(-x*x)+cos(x)+sin(x)+tan(x));
}

int main(int argc, char **argv){
	
	int i, j, k;
	int nRank, nProcs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	MPI_Request req, req2, req3, req4;
	MPI_Status status;

	int n1,n2,jsta,jend;
	int iter,niter;
	double xi,xf,dx;
	double tmr,temp;
	double *ar, *br;
	/* Do not change */
	n1 = 0;
	n2 = 100000000;
	niter = 3;
	/* Do not change */

	int lastR, start, end, remain, jobs;

	lastR = nProcs - 1;
	jobs = n2 / nProcs; //JOB 바꾸삼
	remain = n2 % nProcs; //JOB 바꾸삼
	
	start = nRank * jobs;

	if (nRank == lastR) {
		end = n2;
	}
	else {
		end = start + jobs;
	}

	ar = (double*) malloc(sizeof(double)*n2);
	br = (double*) malloc(sizeof(double)*n2);

	xi = 0.L;
	xf = 1.L;
	dx = (xf-xi)/(double)(n2-n1-1);

	for(i=start;i<end;i++){   // 0~999999데이터 초기화
		br[i] = xi+(double)(i-n1)*dx;
	}
	
	int capyS, capyE;
	capyS = start;
	capyE = end;

	if (nRank == ROOT) {
		start++;
		ar[0] = 0.0;
	}
	if (nRank == lastR) {
		end--;
		ar[end] = 0;
	}

	int endD, startD;
	startD = start - 1;
	endD = end - 1;

	int lastN, nextN;
	lastN = nRank - 1;
	nextN = nRank + 1;
	if (nRank == ROOT)
		lastN = MPI_PROC_NULL;
	if (nRank == lastR)
		nextN = MPI_PROC_NULL;

	for(iter=0;iter<niter;iter++){
		if (nRank == ROOT){
			MPI_Isend(&br[endD], 1, MPI_DOUBLE, nextN, 77, MPI_COMM_WORLD, &req2);
			MPI_Irecv(&br[end], 1, MPI_DOUBLE, nextN, 77, MPI_COMM_WORLD, &req4);
			MPI_Wait(&req2, &status);
			MPI_Wait(&req4, &status);
		}
		else if (nRank == lastR){
			MPI_Isend(&br[start], 1, MPI_DOUBLE, lastN, 77, MPI_COMM_WORLD, &req);
			MPI_Irecv(&br[startD], 1, MPI_DOUBLE, lastN, 77, MPI_COMM_WORLD, &req3);
			MPI_Wait(&req, &status);
			MPI_Wait(&req3, &status);
		}
		else{
			MPI_Isend(&br[start], 1, MPI_DOUBLE, lastN, 77, MPI_COMM_WORLD, &req);
			MPI_Isend(&br[endD], 1, MPI_DOUBLE, nextN, 77, MPI_COMM_WORLD, &req2);
			MPI_Irecv(&br[startD], 1, MPI_DOUBLE, lastN, 77, MPI_COMM_WORLD, &req3);
			MPI_Irecv(&br[end], 1, MPI_DOUBLE, nextN, 77, MPI_COMM_WORLD, &req4);
			MPI_Wait(&req, &status);
			MPI_Wait(&req2, &status);
			MPI_Wait(&req3, &status);
			MPI_Wait(&req4, &status);
		}
		for(j=start;j<end;j++){
			/* Do not change */
			ar[j] = (br[j-1]+br[j+1])/4.L+br[j]/2.L+1.L/genvv(br[j]);
			/* Do not change */
		}
		for(i=capyS;i<capyE;i++){
			/* Do not change */
			br[i] = ar[i];
			/* Do not change */
		}
	}
	temp = 0.L;
	for(j=start;j<end;j++){
		temp += br[j];
	}
	MPI_Reduce(&temp, &tmr, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
	if (nRank == ROOT){
		printf("tmr = %16.7f\n", tmr);
	}
	
	free(ar);
	free(br);
	MPI_Finalize();
	return 0;

}
