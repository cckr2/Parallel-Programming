#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#define ROOT 0

double dot(int N, const double* u, const double* v) {
	double sum=0.0;
	int i;
	for(i=0;i<N;i++) sum+= u[i]*v[i];
	return sum;
}

void PCGM(int N, double** A, double* b, double *res, double abstol) {
	FILE* fp = fopen("C_Serial.txt", "wt");
	

	int i, j, k;
	int nRank, nProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	MPI_Request req;
	MPI_Status status;

	int lastR, start, end, remain, jobs;
	int *displs = (int*)malloc(sizeof(int)*nProcs);
	int *count = (int*)malloc(sizeof(int)*nProcs);

	lastR = nProcs - 1;
	jobs = N / nProcs; //JOB ¹Ù²Ù»ï
	remain = N % nProcs; //JOB ¹Ù²Ù»ï
	start = nRank * jobs;

	if (nRank == lastR) {
		end = start + jobs + remain;
	}
	else {
		end = start + jobs;
	}

	for (i = 0; i < nProcs; i++){
		if (i == lastR) {
			count[i] = jobs + remain;
			displs[i] = jobs*i + remain;
		}
		else {
			count[i] = jobs;
			displs[i] = jobs*i;
		}
	}
	int mycount = count[nRank];
	int maxiters = 100000; 
	int iter;
	double sum, c, d, alpha, beta;
	double tsum, td;

	double *p, *r, *mr, *z,*x;
	double *tp,*tx;
	tp = (double*)malloc(sizeof(double)*mycount);
	x = (double*)malloc(sizeof(double)*N);
	tx = (double*)malloc(sizeof(double)*mycount);
	p = (double*)malloc(sizeof(double)*N);
	r=(double*)malloc(sizeof(double)*N);
	mr=(double*)malloc(sizeof(double)*N);
	z=(double*)malloc(sizeof(double)*N);
	
	//Initial guess
	for (i = 0; i<mycount; i++)
		tx[i] = 1.0;

	for(i=0;i<N;i++)
		p[i] = 1.0;

	for(i=0;i<N;i++) {
		r[i] = b[i] - dot(N, A[i], p);
		mr[i] = r[i]/A[i][i];
	}
	c=dot(N, mr, r);

	for(i=0;i<N;i++) 
		p[i]=mr[i];

	for (iter = 0; iter < maxiters; iter++) {

		tsum = 0.0;
		for (i = start; i < end; i++){
			z[i] = dot(N, A[i], p);
			tsum += p[i] * z[i];
		}
		MPI_Allreduce(&tsum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		/*
		if (nRank == ROOT)
			printf("sum : %.10f\n", sum);
			*/
		MPI_Barrier(MPI_COMM_WORLD);
		alpha = c / sum;
		
		tsum = 0.0;
		td = 0.0;
		j = 0;
		for (i = start; i < end; i++) {
			tx[j] = tx[j] + alpha*p[i];
			r[i] = r[i] - alpha*z[i];
			mr[i] = r[i] / A[i][i];
			td += mr[i] * r[i];
			tsum += r[i] * r[i];
			j++;
		}
		
		MPI_Allreduce(&td, &d, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&tsum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
		/*
		if (nRank == ROOT){
			printf("d : %.10f\n", d);
			printf("sum2 : %.10f\n", sum);
		}
		*/
		if (fabs(d) < abstol) {
			*res = sqrt(fabs(d));
			break;
		}
		if (fabs(sum) < abstol) {
			*res = sqrt(fabs(sum));
			break;
		}
		beta = d / c;

		j = 0;
		for (i = start; i < end; i++){
			tp[j] = mr[i] + beta*p[i];
			j++;
		}
		
		MPI_Allgatherv(tp, mycount, MPI_DOUBLE, p, count, displs, MPI_DOUBLE, MPI_COMM_WORLD);

		c=d;
	}
	MPI_Gatherv(tx, mycount, MPI_DOUBLE, x, count, displs, MPI_DOUBLE,0, MPI_COMM_WORLD);
	if (nRank == ROOT){
		printf("iters = %d, residual norm=%g\n", iter+1, res);
		for (i = 0; i < N; i++) {
			fprintf(fp, "%d		%lf\n", i, x[i]);
		}
	}
	fclose(fp);
	free(x);
	free(tx);

	free(p); free(z); free(r); free(mr);
	
}

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	int nRank, nProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	int i, j;
	int iters;
	int const N = 3200;
	double **A,  *b;
	double res;

	b=(double*)malloc(sizeof(double)*N);
	
	A=(double**)malloc(sizeof(*A)*N);
	for(i=0;i<N;i++) A[i]=(double*)malloc(sizeof(double)*N);

	for(i=0;i<N;i++) {
		b[i]=i+1;
		for(j=0;j<N;j++) {
			A[i][j] = 0.0;
		}
	}
	for(i=0;i<N;i++) {
		A[i][i] = -2.0;
		if(i<N-1) {
			A[i][i+1] = 1.0;
			A[i+1][i] = 1.0;
		}	
	}


	
	PCGM(N, A,  b, &res, 1e-10);
	
	for(i=0;i<N;i++) free(A[i]);
	free(A);
	
	free(b);

	MPI_Finalize();
	return 0;
}
