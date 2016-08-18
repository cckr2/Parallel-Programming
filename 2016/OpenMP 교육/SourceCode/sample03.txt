#include <stdio.h>
#include <omp.h>
int main()
{
    int i = 10, tid;
    omp_set_num_threads(4);
    printf("start tid = X i=%d\n", i);
#pragma omp parallel private(tid) firstprivate(i)
    {
        tid = omp_get_thread_num();
        printf ("tid = %d i = %d\nr", tid, i);
        i = 20;
    }
    printf("end tid = %d i=%d\n", tid, i);
    return 0;
}
