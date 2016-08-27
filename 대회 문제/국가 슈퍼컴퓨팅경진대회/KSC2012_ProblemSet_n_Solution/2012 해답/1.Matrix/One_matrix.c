#include <stdio.h>
#include <sys/time.h>

#define N 512

int timeval_subtract(struct timeval *, struct timeval *, struct timeval *);

int main()
{
    int i, j, k;
    double a[N][N], b[N][N], c[N][N];
    double sum=0.0, flops=0.0;

    struct timeval tv_start, tv_end, tv_diff;

    for(i=0;i<N;i++)
	for(j=0;j<N;j++)
	    c[i][j]=0.0;

    printf("OK\n");

    gettimeofday(&tv_start, NULL);

    for(i=0;i<N;i++) {
	for(j=0;j<N;j++) {
	    a[i][j] = (double)(i+1)/(double)(j+1);
	    b[i][j] = (double)(i+1)/(double)(j+1);
	}
    }

    for(i=0;i<N;i++)
	for(j=0;j<N;j++)
		for(k=0;k<N;k++)
		    c[i][j] += a[i][k]*b[k][j];

    for(i=0;i<N;i++)
	for(j=0;j<N;j++)
	    sum += sqrt(c[i][j]*c[i][j]); // ??????

    gettimeofday(&tv_end, NULL);

    timeval_subtract(&tv_diff, &tv_end, &tv_start);
    printf("Elapsed time : %ld.%06ld ms. sum : %lf\n", tv_diff.tv_sec, tv_diff.tv_usec, sum);

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
