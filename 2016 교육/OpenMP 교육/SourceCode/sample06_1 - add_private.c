#include <stdio.h>
#include <math.h>
#define num_steps 1000000000
#define NUM_THREADS 8

int main(int argc, char *argv[]) {
    double sum[NUM_THREADS], step, x, pi;
    double t1, t2;
    int i, id;

    step=1./(double)num_steps;

#pragma omp parallel private(id,x)
    {
        id = omp_get_thread_num();
        sum[id]=0.0;
        for (i=id; i< num_steps; i=i+NUM_THREADS){
            x = (i+0.5)*step;
            sum[id] += 4.0/(1.0+x*x);
        }
    }
    for(i=0, pi=0.0;i<NUM_THREADS;i++)
        pi += sum[i] * step;

    printf(" numerical pi = %.15f \n", pi);
    printf(" analytical pi = %.15f \n", acos(-1.0));
    printf(" Error = %E \n", fabs(acos(-1.0)-pi));
    return 0;
}
