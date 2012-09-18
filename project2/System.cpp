#include <iostream>
#include "math.h"
#include "time.h"
#include <fstream>
#include "Atom.h"
#include "System.h"

using namespace std;

inline float squaredDistanceBetweenAtoms(Atom *atom1, Atom *atom2) {
	return norm((atom1->r - atom2->r),1);
}

inline vec vectorBetweenAtoms(Atom *atom1, Atom *atom2) {
	return atom1->r - atom2->r;
}

inline vec forceBetweenAtoms(Atom *atom1, Atom *atom2) {
	return zeros<vec>(3,1);
}

System::System(int N, double T, double rho) {
	this->N = N;
	this->T = T;
	this->rho = rho;
	this->initialize();
}

void System::step(double dt) {
	for(int n=0;n<this->N;n++) 
		this->atoms[n]->a.zeros();

	vec F = zeros<vec>(3,1);
	Atom *atom0, *atom1;
	float rsq, invSqrt;

	mat v_temp = zeros<mat>(3,this->N);
	
	for(int i=0;i<this->N-1;i++) {
		for(int j=i+1;j<this->N;j++) {
			atom0 = this->atoms[i];
			atom1 = this->atoms[j];

			rsq = squaredDistanceBetweenAtoms(atom0,atom1);
			invSqrt = FastInvSqrt(rsq);

			F = forceBetweenAtoms(atom0,atom1);
			atom0->a += F/atom0->mass;
			atom1->a -= F/atom1->mass;
		}
	}

	for(int n=0;n<this->N;n++) {
		v_temp.col(n) = this->atoms[n]->v + 0.5*this->atoms[n]->a*dt;
		this->atoms[n]->r += v_temp.col(n)*dt;
	}

	for(int n=0;n<this->N;n++) 
		this->atoms[n]->a.zeros();

	for(int i=0;i<this->N-1;i++) {
		for(int j=i+1;j<this->N;j++) {
			atom0 = this->atoms[i];
			atom1 = this->atoms[j];

			rsq = squaredDistanceBetweenAtoms(atom0,atom1);
			invSqrt = FastInvSqrt(rsq);

			F = forceBetweenAtoms(atom0,atom1);
			atom0->a += F/atom0->mass;
			atom1->a -= F/atom1->mass;
		}
	}

	for(int n=0;n<this->N;n++) {
		this->atoms[n]->v = v_temp.col(n) + 0.5*this->atoms[n]->a*dt;
	}
}


void System::initialize() {
	this->atoms = new Atom*[this->N];
	for(int n=0;n<this->N;n++) {
		this->atoms[n] = new Atom();
	}

	this->initPositions();
	this->initVelocities();
}

void System::initPositions() {
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
}

double System::gasdev() {
	static bool available = false;
	static double gset;
	double fac, rsq, v1, v2;
	if(!available) {
		do {
			v1 = 2.0*rand() / double(RAND_MAX) - 1.0;
			v2 = 2.0*rand() / double(RAND_MAX) - 1.0;
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
	for(int n=0;n<this->N;n++) 
		for(int i=0;i<3;i++) 
			this->atoms[n]->v(i) = this->gasdev();

	vec vCM = zeros<vec>(3,1);
	
	for(int n=0;n<this->N;n++)
		vCM += this->atoms[n]->v;
		
	vCM /= this->N;
	
	for(int n=0;n<this->N;n++) 
		this->atoms[n]->v -= vCM;

	this->rescaleVelocities();
}

void System::rescaleVelocities() {
	double vSqdSum = 0;
	for(int n=0;n<this->N;n++)
		vSqdSum += norm(this->atoms[n]->v,1);
	double lambda = sqrt(3*(this->N-1)*this->T/vSqdSum);
	for(int n=0;n<this->N;n++) 
		this->atoms[n]->v *= lambda;
}

Atom ** System::getAtoms() { return this->atoms; }
int System::getNumberOfAtoms() { return this->N; }
double System::getTemperature() { return this->T; }
double System::getDensity() { return this->rho; }
double System::getLength() { return this->L; }

void System::printVelocitiesToScreen() {
	for(int n=0;n<this->N;n++) {
		printf("r[%d]=(%f,%f,%f)\n",n,this->atoms[n]->v(0),this->atoms[n]->v(1),this->atoms[n]->v(2));
	}
}

void System::printPositionsToFile(ofstream *file) {
	*file << this->N << endl;
	*file << "H atoms are the face atoms, O are the cube atoms" << endl;

	for(int n=0;n<this->N;n++) {
		*file << (this->atoms[n]->type ? "H " : "O ") << this->atoms[n]->r(0) << " " << this->atoms[n]->r(1) << " " << this->atoms[n]->r(2) << endl;
	}
	
}

float FastInvSqrt(float x) {
  float xhalf = 0.5f * x;
  int i = *(int*)&x;         // evil floating point bit level hacking
  i = 0x5f3759df - (i >> 1);  // what the fuck?
  x = *(float*)&i;
  x = x*(1.5f-(xhalf*x*x));
  return x;
}