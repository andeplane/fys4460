#define VERLET_LISTS
// #define RESCALE_VELOCITIES

#include <iostream>
#include "math.h"
#include "time.h"
#include <fstream>
#include "Atom.h"
#include "System.h"
#include "lib.h"
#include "initialConditions.cpp"

using namespace arma;
using namespace std;

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
	this->P = this->rho*this->T;
	double volume = pow(this->L,3.0);

	for(int n=0;n<this->N;n++) {
		this->atoms[n]->F.zeros();
		this->atoms[n]->a.zeros();	
		this->atoms[n]->potential_energy = 0;
	}

	vec dr = zeros<vec>(3,1);
	double dr_2;
	double dr_6;
	double dr_12;

	Atom *atom0;
	Atom *atom1;
	double f,potential_energy;
	
	double dP = 0;
	
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
			dP += f*norm(dr,2);

			potential_energy = 4*(1.0/dr_12 - 1.0/dr_6);

			atom0->a += f*dr;
			atom0->potential_energy += potential_energy;
			atom1->a -= f*dr;
			atom1->potential_energy += potential_energy;
		}
	}
	this->P += 1.0/(3*volume)*dP;

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

void System::printPositionsToFile(ofstream *file) {
	*file << this->N << endl;
	*file << "H atoms are the face atoms, O are the cube atoms" << endl;

	for(int n=0;n<this->N;n++) {
		*file << (this->atoms[n]->type ? "H " : "O ") << this->atoms[n]->r(0) << " " << this->atoms[n]->r(1) << " " << this->atoms[n]->r(2) << endl;
	}
	
}