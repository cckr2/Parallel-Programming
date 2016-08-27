#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include "mpi.h"

#define im 200
#define jm 200
#define im1 201
#define jm1 201
#define ITER_MAX 100000

int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y);

int main(int argc, char* argv[])
{
    int nprocs, nrank;
    int start, end, prev, next;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &nrank);
 
    MPI_Request isend1, irecv1, isend2, irecv2;
    MPI_Status status;

    int i, j, iter, nprt;

    double tolerance, error, t_error;
    double u[im1][jm1], uo[im1][jm1];
    double ub[im1][jm1], ut[im1][jm1];
    double elapsed_time;

    struct timeval tv_start, tv_end, tv_diff;

    double *buff_prev, *buff_next;

    //Set MPI vars
    start = im/nprocs*nrank;
    end = im/nprocs*(nrank+1)-1;
    
    start = (!nrank) ? 1 : start;

    prev = (nrank) ? nrank-1 : -1;
    next = (nrank==nprocs-1) ? -1 : nrank+1;

    //Read Input
    iter=0;
    nprt=500;
    tolerance = 1e-6;
    elapsed_time = 0.0;
    
    for(i=0;i<im1;i++) {
	for(j=0;j<jm1;j++) {
	    u[i][j]=0.0;
	    uo[i][j]=0.0;
	    ub[i][j]=0.0;
	    ut[i][j]=0.0;
	}
    }

    buff_prev = (double *)calloc(jm1,sizeof(double));
    buff_next = (double *)calloc(jm1,sizeof(double));

    //Boundary condition
    for(j=0; j<=jm; j++) {
	u[0][j] = 1.0;
	u[im][j] = 1.0;
    }
    for(i=0; i<=im; i++) {
	u[i][0] = 1.0;
	u[i][jm] = 2.0;
    }

    gettimeofday(&tv_start, NULL);
    
    do {

	for(i=0; i<=im; i++)
	    for(j=0; j<=jm; j++)
		uo[i][j] = u[i][j];
	
	for(j=0; j<=jm; j++) {
	    buff_prev[j] = uo[start][j];
	    buff_next[j] = uo[end][j];
	}

	//Communication
	if(prev != -1) MPI_Isend(buff_prev, jm1, MPI_DOUBLE, prev, 10, MPI_COMM_WORLD, &isend1);
	if(next != -1) MPI_Irecv(u[end+1], jm1, MPI_DOUBLE, next, 10, MPI_COMM_WORLD, &irecv1);
	
	if(next != -1) MPI_Isend(buff_next, jm1, MPI_DOUBLE, next, 11, MPI_COMM_WORLD, &isend2);
	if(prev != -1) MPI_Irecv(u[start-1], jm1, MPI_DOUBLE, prev, 11, MPI_COMM_WORLD, &irecv2);

	if(prev != -1) {
	    MPI_Wait(&isend1, &status);
	    MPI_Wait(&irecv2, &status);
	}
	if(next != -1) {
	    MPI_Wait(&irecv1, &status);
	    MPI_Wait(&isend2, &status);
	}

	for(i=start; i<=end; i++)
	    for(j=1; j<jm; j++)
		u[i][j] = (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1])/4;

	error = 0.0;
	for(i=start; i<=end; i++)
	    for(j=0; j<=jm; j++)
		error += ((u[i][j]-uo[i][j])*(u[i][j]-uo[i][j]));

	MPI_Allreduce(&error, &t_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    } while((iter++ < ITER_MAX) && (t_error > tolerance));


    gettimeofday(&tv_end, NULL);
    
    for(i=start; i<=end; i++)
	for(j=1; j<jm; j++)
	    ub[i][j]=u[i][j];
    
    MPI_Reduce(ub, ut, im1*jm1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
    
    for(j=0; j<=jm; j++) {
	ut[0][j] = 1.0;
	ut[im][j] = 1.0;
    }
    for(i=0; i<=im; i++) {
	ut[i][0] = 1.0;
	ut[i][jm] = 2.0;
    }

    timeval_subtract(&tv_diff, &tv_end, &tv_start);
    
    if(nrank==0) {
	printf("Elapsed time : %ld.%06ld s. iter : %d, error : %e\n", tv_diff.tv_sec, tv_diff.tv_usec,iter,t_error);

	FILE *fp;

	fp = fopen("FDM2D_parallel.txt","w");

	for(i=0;i<im1;i++) {
	    for(j=0;j<jm1;j++) {
		fprintf(fp,"%d %d %.4lf\n",i,j,ut[i][j]);
	    }
	}

	fclose(fp);
    }

    MPI_Finalize();
    return 0;
}

int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y)
{
    int nsec;

    if(x->tv_usec < y->tv_usec) {
	nsec = (y->tv_usec - x->tv_usec)/1000000 + 1;
	y->tv_usec -= 1000000*nsec;
	y->tv_sec += nsec;
    }
    if(x->tv_usec - y->tv_usec > 1000000) {
	nsec = (y->tv_usec - x->tv_usec)/1000000;
	y->tv_usec += 1000000*nsec;
	y->tv_sec-=nsec;
    }
    result->tv_sec = x->tv_sec - y->tv_sec;
    result->tv_usec = x->tv_usec - y->tv_usec;

    return 0;
}
