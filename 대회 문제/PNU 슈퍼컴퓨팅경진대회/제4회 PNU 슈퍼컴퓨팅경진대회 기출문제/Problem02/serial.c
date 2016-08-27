#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double dot(int N, const double* u, const double* v) {
	double sum=0.0;
	int i;
	for(i=0;i<N;i++) sum+= u[i]*v[i];
	return sum;
}

int PCGM(int N, double** A, double* x, double* b, double *res, double abstol) {
	int i,j,k;
	int maxiters = 2;
	int iter;
	double sum, c, d, alpha, beta;

	double *p, *r, *mr, *z;
	p=(double*)malloc(sizeof(double)*N);
	r=(double*)malloc(sizeof(double)*N);
	mr=(double*)malloc(sizeof(double)*N);
	z=(double*)malloc(sizeof(double)*N);
	
	//Initial guess
	for(i=0;i<N;i++) p[i] = 1.0;
	for(i=0;i<N;i++) {
		r[i] = b[i] - dot(N, A[i], p);
		mr[i] = r[i]/A[i][i];
	}
	c=dot(N, mr, r);

	for(i=0;i<N;i++) p[i]=mr[i];

	for(iter=0;iter<maxiters;iter++) {
		for(i=0;i<N;i++) z[i]=dot(N,A[i],p);
/*		if(iter==1) {
		//	printf("sum = %.30f\n", sum);
			for(i=0;i<N;i++) {
				printf("z[%d] = %.30f\n", i, z[i]);
			}
		}*/
		sum=dot(N,p,z);
		printf("sum : %.10f\n", sum);
		alpha=c/sum;
		printf("alpha : %f\n", alpha);
		for(i=0;i<N;i++) {
			x[i]=x[i]+alpha*p[i];
			r[i]=r[i]-alpha*z[i];
		}
		/* Preconditioning Stage */
		for(i=0;i<N;i++) 
			mr[i] = r[i]/A[i][i];
		
		d=dot(N,mr,r);
		sum=dot(N,r,r);
		printf("d : %.10f\n", d);
		printf("sum2 : %.10f\n", sum);
		if(fabs(d)<abstol) {
			*res=sqrt(fabs(d));
			break;
		}
		if(fabs(sum)<abstol) {
			*res=sqrt(fabs(sum));
			break;
		}
		beta=d/c;
		for(i=0;i<N;i++) p[i]=mr[i]+beta*p[i];
		c=d;
	}
	free(p); free(z); free(r); free(mr);
	return iter+1;
}

int main(int argc, char* argv[]) {
	int i, j;
	int iters;
	int const N = 3200;
	double **A, *x, *b;
	double res;

	FILE* fp = fopen("C_Serial.txt", "wt");
	x=(double*)malloc(sizeof(double)*N);
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
	for(i=0;i<N;i++)
		x[i] = 1.0;

	iters = PCGM(N, A, x, b, &res, 1e-10);

	printf("iters = %d, residual norm=%g\n",iters, res);

	for(i=0;i<N;i++) {
		fprintf(fp, "%d		%lf\n",i,x[i]);
	}
	fclose(fp);

	for(i=0;i<N;i++) free(A[i]);
	free(A);
	free(x);
	free(b);
	return 0;
}
