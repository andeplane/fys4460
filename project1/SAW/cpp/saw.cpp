#include <iostream>
#include "math.h"
#include "lib.h"
#include "time.h"
#include <armadillo>
#define SHOWPERCENTAGE
using namespace arma;

inline void func(int a) {
	
}

int main(int argc, char *argv[]) {
	int k_max = 9;
	int N_max = pow(2.0,k_max);
	int walkers = atoi(argv[1]);

	long idum;
	
	idum = -1;
	
	mat A = randu<mat>(1,1);
	A.zeros(2*N_max+1,2*N_max+1);
	
	mat rabs = randu<mat>(k_max,walkers);
	rabs.zeros();

	mat r = randu<mat>(1,1);
	r.zeros(N_max,2);
	
	int u = 0;
	int N;
	int c = 0;
	int cmax = k_max*walkers;
	float maxR = 0;
	
	for(int k=1;k<=k_max;k++) {
		N = pow(2.0,k);

		for(int i=0;i<walkers;i++) {
			
#ifdef SHOWPERCENTAGE
			if(!(++c % (cmax/50))) {
				printf("%d%%...",100*c/cmax);
				fflush(stdout);
			}
#endif
			A.zeros();
			A(N_max,N_max) = 1;

			for(int n=1;n<N_max;n++) {
				u = 4*ran0(&idum);
				
				r(n,0) = r(n-1,0) + (u==0) - (u==1);
				r(n,1) = r(n-1,1) + (u==2) - (u==3);
				
				// printf("(%d,%d)->(%d,%d)\n",(int)r(n-1,0),(int)r(n-1,1),(int)r(n,0),(int)r(n,1));
				
				if(A(r(n,0)+N_max,r(n,1)+N_max) || n>N) {
					// printf("A(%d,%d): %d / n(%d): %d\n",(int)(r(n,0)),(int)(r(n,1)),(int)A(r(n,0)+N_max,r(n,1)+N_max),n,n>N);
					r(n,0) = r(n-1,0);
					r(n,1) = r(n-1,1);

					rabs(k-1,i) = sqrt(r(n,0)*r(n,0)+r(n,1)*r(n,1));
					if(maxR < rabs(k-1,i)) {
						maxR = rabs(k-1,i);
						printf("New maxR: %f\n",maxR);
					}
					break;
				}

				A(r(n,0)+N_max,r(n,1)+N_max) = 1;
			}
		}
	}
	rabs.save("data.dat",raw_ascii);

	return 0;
}