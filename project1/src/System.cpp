// #define RESCALE_VELOCITIES
#define FAST_COLLISIONS

#include <iostream>
#include "math.h"
#include "time.h"
#include <fstream>
#include "Atom.h"
#include "System.h"
#include "InitialConditions.cpp"
#include <inlines.h>

using namespace arma;
using namespace std;

double t = 0;
int steps = 0;

System::System(int number_of_FCC_cells_, double T_, double rho_) {
    number_of_FCC_cells = number_of_FCC_cells_;
    T = T_;
    rho = rho_;
    initialize();
}

double calculate_force_between_atoms(Atom *atom0, Atom *atom1) {
    double dr_2, dr_6, dr_12, f, potential_energy;
    vec dr;

    dr = atom0->distanceToAtom(atom1);
    dr_2 = dot(dr,dr);

    dr_6 = pow(dr_2,3);
    dr_12 = pow(dr_6,2);

    f = 24*(2.0/dr_12-1.0/dr_6)/dr_2;

    potential_energy = 2*(1.0/dr_12 - 1.0/dr_6);

    atom0->a += f*dr;
    atom0->potential_energy += potential_energy;
    atom1->a -= f*dr;
    atom1->potential_energy += potential_energy;

    return f*norm(dr,2);
}

void System::calculateAccelerations() {
    P = rho*T;
    double volume = pow(L,3.0);

    for(int n=0;n<N;n++) {
        atoms[n]->F.zeros();
        atoms[n]->a.zeros();
        atoms[n]->potential_energy = 0;
	}

#ifdef FAST_COLLISIONS
    for(int i=0;i<cells.size();i++)
        cells[i]->reset();

    for(int i=0;i<cells.size();i++) {
        P += 1.0/(3*volume)*cells[i]->calculate_forces(this);
    }
#else
    Atom *atom0, *atom1;
    for(int i=0;i<N;i++) {
        atom0 = atoms[i];
        for(int j=i+1;j<N;j++) {
            atom1 = atoms[j];

            P += calculate_force_between_atoms(atom0,atom1);
        }
    }
#endif
}

void System::step(double dt) {
#ifdef FAST_COLLISIONS
    sort_cells();
#endif

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
        cells[i]->reset_atom_list();
    }

    for(int n=0;n<N;n++) {
        a = atoms[n];

        i = a->r(0)/cell_width;
        j = a->r(1)/cell_width;
        k = a->r(2)/cell_width;

        int cell_index = calculate_cell_index(i,j,k,cells_x,cells_y,cells_z);

        cells[cell_index]->add_atom(a);
    }
}

void System::printPositionsToFile(ofstream *file) {
    *file << N << endl;
	*file << "H atoms are the face atoms, O are the cube atoms" << endl;

    for(int n=0;n<N;n++) {
        *file << (atoms[n]->type ? "H " : "O ") << atoms[n]->r(0) << " " << atoms[n]->r(1) << " " << atoms[n]->r(2) << endl;
	}
	
}
