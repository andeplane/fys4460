// #define VERLET_LISTS
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

void System::calculateAccelerations() {
    P = rho*T;
    double volume = pow(L,3.0);

    for(int n=0;n<N;n++) {
        atoms[n]->F.zeros();
        atoms[n]->a.zeros();
        atoms[n]->potential_energy = 0;
	}

    for(int i=0;i<cells.size();i++) {
        P += 1.0/(3*volume)*cells[i].calculate_forces();
    }
}

void System::step(double dt) {
    sort_cells();
    // time_t t0 = clock();
    // cout << "Time spent on sorting: " << ((double)clock()-t0)/CLOCKS_PER_SEC << endl;

    for(int n=0;n<N;n++) {
        atoms[n]->v += 0.5*atoms[n]->a*dt;
        atoms[n]->addR(atoms[n]->v*dt);
	}

    calculateAccelerations();
    for(int n=0;n<N;n++) {
        atoms[n]->v += 0.5*atoms[n]->a*dt;
	}

    t += dt;
    steps++;

#ifdef RESCALE_VELOCITIES
	if(!(steps % 200)) {
        rescaleVelocities();
	}
#endif
}

void System::sort_cells() {
    int i,j,k;
    Atom *a;
    for(int i=0;i<cells.size();i++) {
        cells[i].reset_atom_list();
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
