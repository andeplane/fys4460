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

double t = 0;
int steps = 0;

System::System(int myid_, Settings *settings_) {
    myid = myid_;
    settings = settings_;
    dt = settings->dt;
    mdio = new MDIO();
    mdio->setup(this);

    rnd = new Random(-(myid+1));
    double b = 1.54478707783;

    Lx = settings->nodes_x*settings->unit_cells_x*b;
    Ly = settings->nodes_y*settings->unit_cells_y*b;
    Lz = settings->nodes_z*settings->unit_cells_z*b;
    if(myid==0) cout << "The system will be of size " << Lx << " " << Ly << " " << Lz << endl;

    this->thread_control = new ThreadControl();
    thread_control->setup(this);
    mdio->save_state_to_file_binary();
    MPI_Finalize();

    exit(0);
    calculateAccelerations();
}

void System::move() {
    for(int i=0;i<thread_control->my_cells.size();i++) {
        Cell *cell = thread_control->my_cells[i];
        Atom *atom = cell->first_atom;
        while(atom != NULL) {
            atom->step(dt);

            int cell_index = thread_control->cell_index_from_atom(atom);
            if(cell_index != atom->cell_index) {
                Cell *new_cell = thread_control->all_cells[cell_index];
                new_cell->new_atoms.push_back(atom); // We will take care of him later
            }
        }
    }

    thread_control->update_cells();
}

void System::calculateAccelerations() {
    /*
    for(unsigned long i=0;i<thread_control->my_cells.size();i++) {
        Cell *cell = thread_control->my_cells[i];
        cell->calculate_forces(this);
    }
    */
}

void System::step() {
    move();
}
