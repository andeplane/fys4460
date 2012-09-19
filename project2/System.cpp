#define VERLET_LISTS
#define RESCALE_VELOCITIES

#include <iostream>
#include "math.h"
#include "time.h"
#include <fstream>
#include "Atom.h"
#include "System.h"
#include "lib.h"

using namespace std;

long idum = -3;
double t = 0;
int steps = 0;

System::System(int N, double T, double rho) {
	this->N = N;
	this->T = T;
	this->rho = rho;
	this->initialize();
}

double maxF = 0;

void System::calculateAccelerations() {
	for(int n=0;n<this->N;n++) {
		this->atoms[n]->F.zeros();
		this->atoms[n]->a.zeros();	
	}

	vec dr = zeros<vec>(3,1);
	double dr_2;
	double dr_6;
	double dr_12;

	Atom *atom0;
	Atom *atom1;
	double f;
	
	for(int i=0;i<this->N-1;i++) {
		atom0 = this->atoms[i];
#ifdef VERLET_LISTS
		for(int j=0;j<atom0->interactingParticles;j++) {
			atom1 = this->atoms[atom0->interactingParticlesList[j]];
#else
		for(int j=i+1;j<this->N;j++) {
			atom1 = this->atoms[j];
#endif
			
			dr = atom0->distanceToAtom(atom1);

			dr_2 = dot(dr,dr);
			
			dr_6 = pow(dr_2,3);
			dr_12 = pow(dr_6,2);
			
			f = 24*(2.0/dr_12-1.0/dr_6)/dr_2;

			atom0->a += f*dr;
			atom1->a -= f*dr;
		}
	}

	vec sumF = zeros<vec>(3,1);
	for(int n=0;n<this->N;n++) {
		sumF += this->atoms[n]->a;
	}
}

void System::step(double dt) {
#ifdef VERLET_LISTS
	this->updateVerletList();
#endif

	t += dt;
	steps++;

	this->calculateAccelerations();
	for(int n=0;n<this->N;n++) {
		this->atoms[n]->v += 0.5*this->atoms[n]->a*dt;
		this->atoms[n]->addR(this->atoms[n]->v*dt + 10*L);
	}

	this->calculateAccelerations();
	for(int n=0;n<this->N;n++) {
		this->atoms[n]->v += 0.5*this->atoms[n]->a*dt;
	}

#ifdef RESCALE_VELOCITIES
	if(!(steps % 200)) {
		this->rescaleVelocities();
	}
#endif
}

void System::updateVerletList() {
	if(steps % 50) return;

	Atom *atom0, *atom1;
	vec dr;
	for(int n=0;n<this->N;n++) 
		this->atoms[n]->interactingParticles = 0;

	for(int i=0;i<this->N-1;i++) {
		atom0 = this->atoms[i];
		for(int j=i+1;j<this->N;j++) {
			atom1 = this->atoms[j];
			dr = atom0->distanceToAtom(atom1);

			if(norm(dr,2) < 3.2) // These will be in contact later
				atom0->interactingParticlesList[atom0->interactingParticles++] = j;
		}
	}
}

void System::initialize() {
	this->atoms = new Atom*[this->N];
	for(int n=0;n<this->N;n++) {
		this->atoms[n] = new Atom(this);
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

void System::printPositionsToFile(ofstream *file) {
	*file << this->N << endl;
	*file << "H atoms are the face atoms, O are the cube atoms" << endl;

	for(int n=0;n<this->N;n++) {
		*file << (this->atoms[n]->type ? "H " : "O ") << this->atoms[n]->r(0) << " " << this->atoms[n]->r(1) << " " << this->atoms[n]->r(2) << endl;
	}
	
}