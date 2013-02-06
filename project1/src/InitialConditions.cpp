#include <System.h>
// #define UNIFORMVELOCITY

double vmax = 8;


void System::initialize() {
    rnd = new Random(-1);

    L = pow(N/rho,1.0/3);
    cells_x = L/3;
    cells_y = L/3;
    cells_z = L/3;

	printf("Initializing system with properties:\n");
    printf("T=%.2f\n",T);
    printf("rho=%.2f\n",rho);
    printf("L=%.2f\n",L);
    printf("Creating %d atoms...",N);

    atoms = new Atom*[N];
    for(int n=0;n<N;n++) {
        atoms[n] = new Atom(this);
	}
	printf("done\n");

    initPositions();
    initVelocities();
    init_cells();
}

void System::init_cells() {
    for(int k=0;k<cells_z;k++) {
        for(int j=0;j<cells_y;j++) {
            for(int i=0;i<cells_x;i++) {
                Cell c;
                cells.push_back(c);
            }
        }
    }
}

void System::initPositions() {
	printf("Initializing FCC lattice...");

	int M=1;
    while(4*M*M*M < N) ++M;
    double a = L/M;

    double xCell[4] = {0, 0.5, 0.5, 0};
    double yCell[4] = {0, 0.5, 0, 0.5};
    double zCell[4] = {0, 0, 0.5, 0.5};
	
	int n = 0;
	for(int x = 0; x < M; x++) {
		for(int y = 0; y < M; y++) {
			for(int z = 0; z < M; z++) {
				for(int k = 0; k < 4; k++) {
					if(n<N) {
						// Set positions and type
                        atoms[n]->r(0) = (x+xCell[k]) * a;
                        atoms[n]->r(1) = (y+yCell[k]) * a;
                        atoms[n]->r(2) = (z+zCell[k]) * a;
                        atoms[n]->r_initial = atoms[n]->r;
						
                        atoms[n]->type = k>0; // For visualization
						++n;
					}
				}
			}
		}
	}
	printf("done\n");
}

double System::gasdev() {
	static bool available = false;
	static double gset;
	double fac, rsq, v1, v2;
	if(!available) {
        do {
            v1 = 2.0*rnd->nextDouble() - 1.0;
            v2 = 2.0*rnd->nextDouble() - 1.0;
			rsq = v1*v1+v2*v2;
        } while(rsq>= 1.0 || rsq == 0.0);

		fac = sqrt(-2.0*log(rsq)/rsq);
		gset = v1*fac;
		available = true;
		return v2*fac;
	}
	else {
		available = false;
		return gset;
	}
}

void System::initVelocities() {
#ifdef UNIFORMVELOCITY
	printf("Creating uniform velocities...");
	double realVmax = vmax/sqrt(3);

    for(int n=0;n<N;n++) {
		for(int k=0;k<3;k++) 
            atoms[n]->v(k) = 2*ran0(&idum)-1;
        atoms[n]->v *= realVmax;
	}
	printf("done\n");
#else

	printf("Creating maxwellian velocities...");

    for(int n=0; n<N; n++ ) {
        atoms[n]->v(0) = rnd->nextGauss()*sqrt(2.0/2*T);
        atoms[n]->v(1) = rnd->nextGauss()*sqrt(2.0/2*T);
        atoms[n]->v(2) = rnd->nextGauss()*sqrt(2.0/2*T);
    }

    /*
    for(int n=0;n<N;n++)
		for(int i=0;i<3;i++) 
            atoms[n]->v(i) = gasdev();
    */

	vec vCM = zeros<vec>(3,1);
	
    for(int n=0;n<N;n++)
        vCM += atoms[n]->v;
		
    vCM /= N;
	
    for(int n=0;n<N;n++)
        atoms[n]->v -= vCM;
	printf("done\n");
    rescaleVelocities();
#endif
}

void System::rescaleVelocities() {
	printf("Rescaling velocities...");
	double vSqdSum = 0;
    for(int n=0;n<N;n++)
        vSqdSum += norm(atoms[n]->v,1);
    double lambda = sqrt(3*(N-1)*T/vSqdSum);
    for(int n=0;n<N;n++)
        atoms[n]->v *= lambda;
	printf("done\n");
}
