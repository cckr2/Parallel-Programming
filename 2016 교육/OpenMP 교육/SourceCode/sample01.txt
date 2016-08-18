#include <stdio.h>
#include <omp.h>
int main(){
#pragma omp parallel
    {
        printf ("Hello World %d\n", omp_get_thread_num());
    }
    printf("\n");
    omp_set_num_threads(4);
#pragma omp parallel
    {
        printf ("Hello World %d\n", omp_get_thread_num());
    }
    printf("\n");
#pragma omp parallel num_threads(2)
    {
        printf ("Hello World %d\n", omp_get_thread_num());
    }
    return 0;
}
