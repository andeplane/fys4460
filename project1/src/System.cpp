#define VERLET_LISTS
// #define RESCALE_VELOCITIES

#include <iostream>
#include "math.h"
#include "time.h"
#include <fstream>
#include "Atom.h"
#include "System.h"
#include "initialConditions.cpp"
#include <inlines.h>

using namespace arma;
using namespace std;

double t = 0;
int steps = 0;

System::System(int N_, double T_, double rho_) {
    N = N_;
    T = T_;
    rho = rho_;
    initialize();
}

double maxF = 0;

void System::calculateAccelerations() {
    P = rho*T;
    double volume = pow(L,3.0);

    for(int n=0;n<N;n++) {
        atoms[n]->F.zeros();
        atoms[n]->a.zeros();
        atoms[n]->potential_energy = 0;
	}

	vec dr = zeros<vec>(3,1);
	double dr_2;
	double dr_6;
	double dr_12;

	Atom *atom0;
	Atom *atom1;
	double f,potential_energy;
	
	double dP = 0;
	
    for(int i=0;i<N-1;i++) {
        atom0 = atoms[i];

#ifdef VERLET_LISTS
		for(int j=0;j<atom0->interactingParticles;j++) {
            atom1 = atoms[atom0->interactingParticlesList[j]];
#else
        for(int j=i+1;j<N;j++) {
            atom1 = atoms[j];
#endif
			
			dr = atom0->distanceToAtom(atom1);

			dr_2 = dot(dr,dr);
			
			dr_6 = pow(dr_2,3);
			dr_12 = pow(dr_6,2);
			
			f = 24*(2.0/dr_12-1.0/dr_6)/dr_2;
			dP += f*norm(dr,2);

			potential_energy = 2*(1.0/dr_12 - 1.0/dr_6);

			atom0->a += f*dr;
			atom0->potential_energy += potential_energy;
			atom1->a -= f*dr;
			atom1->potential_energy += potential_energy;
		}
	}
    P += 1.0/(3*volume)*dP;

	vec sumF = zeros<vec>(3,1);
    for(int n=0;n<N;n++) {
        sumF += atoms[n]->a;
	}
}

void System::step(double dt) {
#ifdef VERLET_LISTS
    updateVerletList();
#endif

	t += dt;
	steps++;

    calculateAccelerations();
    for(int n=0;n<N;n++) {
        atoms[n]->v += 0.5*atoms[n]->a*dt;
        atoms[n]->addR(atoms[n]->v*dt);
	}

    calculateAccelerations();
    for(int n=0;n<N;n++) {
        atoms[n]->v += 0.5*atoms[n]->a*dt;
	}

#ifdef RESCALE_VELOCITIES
	if(!(steps % 200)) {
        rescaleVelocities();
	}
#endif
}

void System::updateVerletList() {
	if(steps % 50) return;

	Atom *atom0, *atom1;
	vec dr;
    for(int n=0;n<N;n++)
        atoms[n]->interactingParticles = 0;

    for(int i=0;i<N-1;i++) {
        atom0 = atoms[i];
        for(int j=i+1;j<N;j++) {
            atom1 = atoms[j];
			dr = atom0->distanceToAtom(atom1);

			if(norm(dr,2) < 3.2) // These will be in contact later
				atom0->interactingParticlesList[atom0->interactingParticles++] = j;
		}
	}
}

void System::sort_cells() {
    int i,j,k;
    Atom *a;
    for(int i=0;i<cells.size();i++) {
        cells[i].reset();
    }

    for(int n=0;n<N;n++) {
        a = atoms[n];
        i = a->r(0)/L;
        j = a->r(1)/L;
        k = a->r(2)/L;
        cells[calculate_cell_index(i,j,k,cells_x,cells_y,cells_z)].add_atom(a);
    }
}

void System::printPositionsToFile(ofstream *file) {
    *file << N << endl;
	*file << "H atoms are the face atoms, O are the cube atoms" << endl;

    for(int n=0;n<N;n++) {
        *file << (atoms[n]->type ? "H " : "O ") << atoms[n]->r(0) << " " << atoms[n]->r(1) << " " << atoms[n]->r(2) << endl;
	}
	
}
