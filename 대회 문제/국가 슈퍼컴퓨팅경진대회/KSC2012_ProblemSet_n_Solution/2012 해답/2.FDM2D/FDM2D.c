#include <stdio.h>
#include <sys/time.h>
#include <math.h>

#define im 200
#define jm 200
#define im1 201
#define jm1 201
#define ITER_MAX 100000

int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y);

int main()
{
    int i, j, iter, nprt;

    double tolerance, error;
    double u[im1][jm1], uo[im1][jm1];
    double elapsed_time;

    struct timeval tv_start, tv_end, tv_diff;

    //Read Input
    iter=0;
    nprt=500;
    tolerance = 1e-6;
    elapsed_time = 0.0;

    //Initialize

    for(i=0; i<=im1; i++) {
	for(j=0; j<=jm1; j++) {
	    u[i][j] = 0.0;
	    uo[i][j] = 0.0;
	}
    }
    
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
	
	for(i=1; i<im; i++)
	    for(j=1; j<jm; j++)
		u[i][j] = (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1])/4;

	error = 0.0;
	for(i=0; i<=im; i++)
	    for(j=0; j<=jm; j++)
		error += ((u[i][j]-uo[i][j])*(u[i][j]-uo[i][j]));

	//printf("%d %.16e\n", iter, error);
    } while((iter++ < ITER_MAX) && (error > tolerance));

    gettimeofday(&tv_end, NULL);
    timeval_subtract(&tv_diff, &tv_end, &tv_start);
    printf("Elapsed time : %ld.%06ld ms. iter : %d, error : %e\n", tv_diff.tv_sec, tv_diff.tv_usec,iter,error);
    
    FILE *fp;

    fp = fopen("FDM2D.txt","w");

    for(i=0;i<im1;i++) {
	for(j=0;j<jm1;j++) {
	    fprintf(fp,"%d %d %.4lf\n",i,j,u[i][j]);
	}
    }

    fclose(fp);

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
