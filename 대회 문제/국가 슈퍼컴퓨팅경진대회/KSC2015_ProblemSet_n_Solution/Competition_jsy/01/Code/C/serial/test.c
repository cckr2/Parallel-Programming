#include <stdio.h>
int main(){
		int iter, niter, j;
			niter = 3;
				int jsta = 0;
					int jend = 100000000 - 1;

						for (iter = 0; iter < niter; iter++){
									for (j = jsta; j < jend; j++){
													/* Do not change */
													printf("ar : %d // br : %d %d %d %d\n",j,j - 1,j + 1,j,j);
															}
										}
}
