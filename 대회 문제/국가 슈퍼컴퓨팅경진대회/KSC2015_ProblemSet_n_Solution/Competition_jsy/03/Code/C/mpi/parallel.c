#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<sys/time.h>
#include<mpi.h>
float ran2(long *);
long iseed=-9;

struct timeval tv;
float gettime(){
	static int startflag = 1;
	static double tsecs0, tsecs1;
	if(startflag) {
		(void ) gettimeofday(&tv, NULL);
		tsecs0 = tv.tv_sec + tv.tv_usec*1.0E-6;
		startflag = 0;
	}
	(void) gettimeofday(&tv, NULL);
	tsecs1 = tv.tv_sec + tv.tv_usec*1.0e-6;
	return (float) (tsecs1 - tsecs0);
}
void domaindecomp(float **ibase, size_t *mmen, float valmin, float valmax, MPI_Comm Comm){
	int niter = 32,temp=0;
	int i,j,nRank;
	MPI_Comm_rank(MPI_COMM_WORLD,&nRank);
	float step = (1.L-0.L)/(float)niter;
	long  number[niter], receive[niter];
	long result;
	printf("%d", mmen);	

	for(i=0;i<mmen;i++){
		j = ibase[i]/step;
		number[j]++;
		temp++;
	}
	for(i=0;i<32;i++){
		if(nRank==i)
		{
			for(j=0;j<32;i++){
				if(nRank==i)
					continue
				else{
					MPI_Irecv(&receive[j], 1, MPI_LONG, j, 20, MPI_COMM_WORLD, &req2);
					MPI_Wait(&req2, &stat);
				}
			}
		}
		else{
			MPI_Isend(&number[i], 1, MPI_LONG, i, 10, MPI_COMM_WORLD, &req1);	
			MPI_Wait(&req1, &stat);
		}
	}
	receive[nRank]=number[nRank];
	for (i=0;i<32;i++){
		result  = result + receive[nRank];
	}
	mmnm = result;
}


int main(int argc, char **argv){
	long i, j;
	long np,niter,maxnp=100000;
	float *val, **lval,lvalmin,lvalmax;
	niter = 32;
	int myid, nid;
	float time1, time2;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	if(nid !=2 && nid !=4 && nid !=8 && nid !=16&&nid!=32){
		fprintf(stderr,"Procs Error\n");
		MPI_Finalize();
	}
	time1 = gettime();
	val = (float*)malloc(sizeof(float)*maxnp);
	float step =(1.-0)/(float)nid;
	iseed = -1*(myid+1284L);
	for(j=0;j<maxnp;j++){
		val[i] = ran2(&iseed);
	}
	np = maxnp;
	float valmin = 0.L; 
	float valmax = 1.L;
	domaindecomp(&val, &np, valmin, valmax, MPI_COMM_WORLD);
	printf("p%d has %ld members\n", myid, np);
	MPI_Barrier(MPI_COMM_WORLD);
	time2 = gettime();
	if (myid == 0){
		printf("Wallclock time = %g second\n", (time2 - time1));
	}
	lvalmin = 1. / nid*(float)(myid);
	lvalmax = 1. / nid*(float)(myid + 1);
	Check(val, np, lvalmin, lvalmax);
	MPI_Finalize();
}
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 4ndecomp(&val, &np, valmin, valmax, MPI_COMM_WORLD);
0692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) {
				iv[j] = *idum;
			}
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
