#pragma once
// #define UNIFORMVELOCITY

long idum = -3;
double vmax = 8;


void System::initialize() {
	printf("Creating %d atoms...",this->N);

	this->atoms = new Atom*[this->N];
	for(int n=0;n<this->N;n++) {
		this->atoms[n] = new Atom(this);
	}
	printf("done\n",this->N);

	this->initPositions();
	this->initVelocities();
}

void System::initPositions() {
	printf("Initiating FCC lattice...");
	this->L = pow(this->N/this->rho,1.0/3);

	int M=1;
	while(4*M*M*M < this->N) ++M;
	double a = this->L/M;
	
	double xCell[4] = {0.25, 0.75, 0.75, 0.25};
	double yCell[4] = {0.25, 0.75, 0.25, 0.75};
	double zCell[4] = {0.25, 0.25, 0.75, 0.75};
	
	int n = 0;
	for(int x = 0; x < M; x++) {
		for(int y = 0; y < M; y++) {
			for(int z = 0; z < M; z++) {
				for(int k = 0; k < 4; k++) {
					if(n<N) {
						// Set positions and type
						this->atoms[n]->r(0) = (x+xCell[k]) * a;
						this->atoms[n]->r(1) = (y+yCell[k]) * a;
						this->atoms[n]->r(2) = (z+zCell[k]) * a;
						this->atoms[n]->type = k>0; // For visualization

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
			v1 = 2.0*ran0(&idum) - 1.0;
			v2 = 2.0*ran0(&idum) - 1.0;
			rsq = v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);

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

	for(int n=0;n<this->N;n++) {
		for(int k=0;k<3;k++) 
			this->atoms[n]->v(k) = 2*ran0(&idum)-1;
		this->atoms[n]->v *= realVmax;
	}
	printf("done\n");
#else
	printf("Creating maxwellian velocities...");
	for(int n=0;n<this->N;n++) 
		for(int i=0;i<3;i++) 
			this->atoms[n]->v(i) = this->gasdev();

	vec vCM = zeros<vec>(3,1);
	
	for(int n=0;n<this->N;n++)
		vCM += this->atoms[n]->v;
		
	vCM /= this->N;
	
	for(int n=0;n<this->N;n++) 
		this->atoms[n]->v -= vCM;
	printf("done\n");
	this->rescaleVelocities();
#endif
}

void System::rescaleVelocities() {
	printf("Rescaling velocities...");
	double vSqdSum = 0;
	for(int n=0;n<this->N;n++)
		vSqdSum += norm(this->atoms[n]->v,1);
	double lambda = sqrt(3*(this->N-1)*this->T/vSqdSum);
	for(int n=0;n<this->N;n++) 
		this->atoms[n]->v *= lambda;
	printf("done\n");
}
