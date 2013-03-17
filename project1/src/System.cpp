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

    Lx = settings->nodes_x*settings->unit_cells_x*settings->FCC_b;
    Ly = settings->nodes_y*settings->unit_cells_y*settings->FCC_b;
    Lz = settings->nodes_z*settings->unit_cells_z*settings->FCC_b;
    this->thread_control = new ThreadControl();
    thread_control->setup(this);
    calculate_accelerations();
    half_kick();

    if(myid==0) cout << "System size: " << Lx << " " << Ly << " " << Lz << endl;
    if(myid==0) cout << "Atoms: " << thread_control->num_atoms << endl;
}
double comm_time = 0;

void System::move() {
    for(unsigned int i=0;i<thread_control->my_cells.size();i++) {
        Cell *cell = thread_control->my_cells[i];
        Atom *atom = cell->first_atom;
        while(atom != NULL) {
            atom->step(dt);

            int cell_index = thread_control->cell_index_from_atom(atom);
            if(cell_index != atom->cell_index) {
                if(cell_index < 0) {
                    cout << "Atom " << atom->index << " is weird. R: " << atom->r[0] << " " << atom->r[1] << " " << atom->r[2] << endl;
                }
                Cell *new_cell = thread_control->all_cells[cell_index];
                new_cell->new_atoms.push_back(atom); // We will take care of him later
            }

            atom = atom->next;
        }
    }
}

void System::half_kick() {
    for(unsigned int n=0;n<thread_control->all_atoms.size();n++) {
        Atom *atom = thread_control->all_atoms[n];
        atom->v[0] += atom->a[0]*dt*0.5;
        atom->v[1] += atom->a[1]*dt*0.5;
        atom->v[2] += atom->a[2]*dt*0.5;
    }
}

void System::kick() {
    for(unsigned int n=0;n<thread_control->all_atoms.size();n++) {
        Atom *atom = thread_control->all_atoms[n];
        atom->v[0] += atom->a[0]*dt;
        atom->v[1] += atom->a[1]*dt;
        atom->v[2] += atom->a[2]*dt;
    }
}

void System::calculate_accelerations() {
    thread_control->reset_forces();

    for(unsigned long i=0;i<thread_control->my_cells.size();i++) {
        Cell *cell = thread_control->my_cells[i];
        cell->calculate_forces(this);
    }
}

void System::step() {
    move();

    double t0 = MPI_Wtime();
    thread_control->update_cells_local();
    thread_control->update_cells_mpi();
    thread_control->update_ghost_cells();
    double t1 = MPI_Wtime();
    comm_time += t1-t0;

    calculate_accelerations();
    // kick();


    steps++;

    if(myid==0 && !(steps % 100)){
        cout << steps <<  endl;
        cout << comm_time << endl;
    }
}
