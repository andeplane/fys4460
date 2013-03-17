#include <iostream>
#include "math.h"
#include "time.h"
#include <fstream>
#include "Atom.h"
#include "System.h"
#include "InitialConditions.cpp"
#include <inlines.h>
#include <settings.h>
#include <mpi.h>
#include <mdio.h>

using namespace arma;
using namespace std;

System::System(int myid_, Settings *settings_) {
    myid = myid_;
    settings = settings_;
    dt = settings->dt;
    mdio = new MDIO();
    mdio->setup(this);
    steps = 0;

    rnd = new Random(-(myid+1));
    double b = 1.54478707783;

    Lx = settings->nodes_x*settings->unit_cells_x*b;
    Ly = settings->nodes_y*settings->unit_cells_y*b;
    Lz = settings->nodes_z*settings->unit_cells_z*b;
    if(myid==0) cout << "The system will be of size " << Lx << " " << Ly << " " << Lz << endl;

    this->thread_control = new ThreadControl();
    thread_control->setup(this);
}

void System::move() {
    double E = 0;
    for(unsigned int i=0;i<thread_control->my_cells.size();i++) {
        Cell *cell = thread_control->my_cells[i];
        Atom *atom = cell->first_atom;
        while(atom != NULL) {
            atom->v[0] += atom->a[0]*dt;
            atom->v[1] += atom->a[1]*dt;
            atom->v[2] += atom->a[2]*dt;
            E += atom->v[2]*atom->v[2] + atom->v[1]*atom->v[1] + atom->v[0]*atom->v[0];

            atom->step(dt);

            int cell_index = thread_control->cell_index_from_atom(atom);
            if(cell_index != atom->cell_index) {
                Cell *new_cell = thread_control->all_cells[cell_index];
                new_cell->new_atoms.push_back(atom); // We will take care of him later
            }

            atom = atom->next;
        }
    }

    thread_control->update_cells_local();
    thread_control->update_cells_mpi();
    // cout << E << endl;
}

long particle_pairs_system = 0;

void calculate_force_between_atoms(Atom *atom0, Atom *atom1, System *system) {
    double dr_2, dr_6, dr_12, f, potential_energy, dr_12_inv, dr_6_inv;
    double dx,dy,dz;
    particle_pairs_system += atom0->index + atom1->index;

    dx = atom0->r[0] - atom1->r[0];
    dy = atom0->r[1] - atom1->r[1];
    dz = atom0->r[2] - atom1->r[2];

    if(dx>system->Lx/2)  dx -= system->Lx;
    if(dx<-system->Lx/2) dx += system->Lx;

    if(dy>system->Ly/2)  dy -= system->Ly;
    if(dy<-system->Ly/2) dy += system->Ly;

    if(dz>system->Lz/2)  dz -= system->Lz;
    if(dz<-system->Lz/2) dz += system->Lz;

    dr_2 = dx*dx + dy*dy + dz*dz;
    // dr_2 = max(0.8,dr_2);

    dr_6 = pow(dr_2,3);
    dr_12 = pow(dr_6,2);
    dr_12_inv = 1.0/dr_12;
    dr_6_inv = 1.0/dr_6;

    f = 12*(2.0*dr_12_inv-dr_6_inv)/dr_2;

    potential_energy = 4*(dr_12_inv - dr_6_inv);

    atom0->a[0] += dx*f;
    atom0->a[1] += dy*f;
    atom0->a[2] += dz*f;

    atom1->a[0] -= dx*f;
    atom1->a[1] -= dy*f;
    atom1->a[2] -= dz*f;
    // particle_pairs_system++;
}

void System::calculate_accelerations() {
    thread_control->reset_forces();

    for(unsigned long i=0;i<thread_control->my_cells.size();i++) {
        Cell *cell = thread_control->my_cells[i];
        cell->calculate_forces(this);
    }
    return;

    for(int i=0;i<thread_control->all_atoms.size();i++) {
        Atom *atom0 = thread_control->all_atoms[i];
        for(int j=i+1;j<thread_control->all_atoms.size();j++) {
            if(i==j) continue;
            Atom *atom1 = thread_control->all_atoms[j];
            calculate_force_between_atoms(atom0,atom1,this);
        }
    }

    cout << particle_pairs_system << endl;
}

void System::step() {
    calculate_accelerations();
    move();
    steps++;
    cout << steps << endl;
}
