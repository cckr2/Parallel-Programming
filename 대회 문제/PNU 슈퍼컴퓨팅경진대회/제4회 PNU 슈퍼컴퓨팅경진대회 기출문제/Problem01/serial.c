#include <stdio.h>
#include<math.h>

#define L 1.0
#define LAMBDA 0.0003
#define PI 3.141592

int main(int argc, char* argv[]){
	long int count = 0, num_modes;
	double x, y, z, r, delta, sol;

	delta = LAMBDA;
	long int test=0;
	for (x = 0.0; x <=L ; x = x + delta){
		printf("%f\n",x);
		/*
		for (y = 0.0; y <= L; y = y + delta){
			for (z = 0.0; z <= L; z = z + delta){
				r = sqrt(x*x + y*y + z*z);
				if (r <= L)
					count++;
			}
		}
		*/
		test++;
	}

	printf("%ld\n",test);
	num_modes = 16 * count;
	sol = 8.0*PI*L *L*L / 3.0 / LAMBDA / LAMBDA / LAMBDA;

	printf("Number of modes = %ld (from mode counting)\n", num_modes);
	printf("Number of modes = %f (from solution)\n", sol);
	
	return 0;
}
