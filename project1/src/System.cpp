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

System::System(int rank_, int nodes_, double dt, int number_of_FCC_cells_, double T_) {
    number_of_FCC_cells = number_of_FCC_cells_;
    T = T_;
    rank = rank_;
    nodes = nodes_;

    initialize(dt);
}

void System::update_velocity_and_move(const double &dt) {
    for(int n=0;n<all_atoms.size();n++) {
        all_atoms[n]->v += 0.5*all_atoms[n]->a*dt;
        all_atoms[n]->addR(all_atoms[n]->v*dt);
    }
}

void System::calculateAccelerations() {
#ifdef MPI_ENABLED
    P = 0;
    if(rank==0) {
        // Reset all atom accelerations
        for(int n=0;n<all_atoms.size();n++) {
            all_atoms[n]->a.zeros();
            all_atoms[n]->potential_energy = 0;
        }
        send_particles_to_slaves();
    }
    else {
        receive_particles_from_master();
    }

    ThreadNode &node = thread_control->nodes[rank];

    for(set<int>::iterator it=node.connected_cells.begin(); it!= node.connected_cells.end();it++) {
        int cell_index = *it;
        cells[cell_index]->reset();
    }

    for(set<int>::iterator it=node.owned_cells.begin(); it!= node.owned_cells.end();it++) {
        int cell_index = *it;
        cells[cell_index]->calculate_forces(this);
    }

    P /= 3*V;

    if(rank == 0) {
        receive_particles_back_from_slaves();
    } else {
        send_particles_back_to_master();
    }
#else
    P = 0;
    for(int n=0;n<N;n++) {
        atoms[n]->a.zeros();
        atoms[n]->potential_energy = 0;
    }

    for(int c=0;c<cells.size();c++) {
        cells[c]->reset();
    }

    for(int c=0;c<cells.size();c++) {
        cells[c]->calculate_forces(this);
    }

    P /= 3*V;
#endif
}

void System::step(double dt) {
    if(rank == 0) {
        update_velocity_and_move(dt);
    }

    calculateAccelerations();

    t += dt;
    steps++;

    if(rank == 0) {
        cout << "Calculating timestep " << steps << endl;
    }
}

void System::sort_cells(bool calculate_all_atoms) {
    int i,j,k;
    Atom *a;
    for(int i=0;i<cells.size();i++) {
        cells[i]->reset_atom_list();
    }

    vector<Atom*> &atom_list = atoms;
    if(calculate_all_atoms) {
        atom_list = all_atoms;
    }

    for(int n=0;n<N;n++) {
        a = atom_list[n];

        i = a->r(0)/cell_width;
        j = a->r(1)/cell_width;
        k = a->r(2)/cell_width;

        int cell_index = calculate_cell_index(i,j,k,cells_x,cells_y,cells_z);

#ifdef DEBUG_ME
        if(cell_index < 0 || cell_index >= cells.size()) {
            cout << "Problem with cell index for particle " << n << endl;
            cout << "Particle is at " << a->r << endl;
            cout << "which gives (i,j,k)=" << i << " " << j << " " << k << endl;
            cout << "Btw, we should have " << N << " atoms." << endl;
        }
#endif

        cells[cell_index]->add_atom(a);
    }
}

void System::printPositionsToFile(ofstream *file) {
    *file << all_atoms.size() << endl;
	*file << "H atoms are the face atoms, O are the cube atoms" << endl;

    for(int n=0;n<all_atoms.size();n++) {
        *file << (all_atoms[n]->type ? "H " : "O ") << all_atoms[n]->r(0) << " " << all_atoms[n]->r(1) << " " << all_atoms[n]->r(2) << endl;
	}
}

#ifdef MPI_ENABLED
void System::send_particles_to_slaves() {
    sort_cells(true);

    double *data = new double[4*N];

    Cell *cell;

    // Add my own atoms to my atoms-list
    int particles = 0;
    N = 0;
    Atom *atom;

    ThreadNode &node = thread_control->nodes[0];
    for(set<int>::iterator it=node.connected_cells.begin(); it!= node.connected_cells.end();it++) {
        int cell_index = *it;
        cell = cells[cell_index];
        for(int n=0;n<cell->atoms.size();n++) {
            atom = atoms[particles];

            if(particles>=atoms.size()) {
                atoms.push_back(atom);
                N = atoms.size();
            }
            else {
                atom = atoms[particles];
            }

            atom->r(0) = cell->atoms[n]->r(0);
            atom->r(1) = cell->atoms[n]->r(1);
            atom->r(2) = cell->atoms[n]->r(2);
            atom->index = cell->atoms[n]->index;
        }
    }

    N = particles;
    particles = 0;

    for(int node_id=1;node_id<nodes;node_id++) {
        ThreadNode &node = thread_control->nodes[node_id];
        for(set<int>::iterator it=node.connected_cells.begin(); it!= node.connected_cells.end();it++) {
            int cell_index = *it;
            cell = cells[cell_index];
            for(int n=0;n<cell->atoms.size();n++) {
                data[particles*4 + 0] = cell->atoms[n]->r(0);
                data[particles*4 + 1] = cell->atoms[n]->r(1);
                data[particles*4 + 2] = cell->atoms[n]->r(2);
                data[particles*4 + 3] = cell->atoms[n]->index;

                particles++;
            }
        }

        MPI_Send(&particles,1,MPI_INT,node_id,100,MPI_COMM_WORLD);
        MPI_Send(data,4*particles,MPI_DOUBLE,node_id,100,MPI_COMM_WORLD);
        particles = 0;
    }

    delete data;
}

void System::receive_particles_from_master() {
    /*
    for(int n=0;n<N;n++) {
        delete atoms[n];
    }

    atoms.clear();
    */

    MPI_Status status;

    int particles = 0;
    MPI_Recv(&particles,1,MPI_INT,0,100,MPI_COMM_WORLD,&status);

    double *data = new double[4*particles];

    MPI_Recv(data,4*particles,MPI_DOUBLE,0,100,MPI_COMM_WORLD,&status);
    Atom *atom;
    for(int n=0;n<particles;n++) {
        if(n>=atoms.size()) {
            atom = new Atom(this);
            atoms.push_back(atom);
            N = atoms.size();
        }
        else {
            atom = atoms[n];
        }

        atom->r(0) = data[4*n + 0];
        atom->r(1) = data[4*n + 1];
        atom->r(2) = data[4*n + 2];

        atom->index = round(data[4*n + 3]);
    }

    N = particles;

    sort_cells();

    delete data;
}

void System::receive_particles_back_from_slaves() {
    MPI_Status status;

    for(int i=1;i<nodes;i++) {
        int particles = 0;
        MPI_Recv(&particles,1,MPI_INT,i,100,MPI_COMM_WORLD,&status);

        double *data = new double[5*particles+1];

        MPI_Recv(data,5*particles+1,MPI_DOUBLE,i,100,MPI_COMM_WORLD,&status);
        for(int j=0;j<particles;j++) {
            atoms[ round(data[5*j+4]) ]->a(0) += data[5*j+0];
            atoms[ round(data[5*j+4]) ]->a(1) += data[5*j+1];
            atoms[ round(data[5*j+4]) ]->a(2) += data[5*j+2];
            atoms[ round(data[5*j+4]) ]->potential_energy += data[5*j+3];
        }

        P += data[5*particles];

        delete data;
    }

    P += N/V*T;
}

void System::send_particles_back_to_master() {
    double *data = new double[5*N+1];

    for(int n=0;n<N;n++) {
        data[5*n+0] = atoms[n]->a(0);
        data[5*n+1] = atoms[n]->a(1);
        data[5*n+2] = atoms[n]->a(2);
        data[5*n+3] = atoms[n]->potential_energy;
        data[5*n+4] = atoms[n]->index;
    }

    data[5*N] = P;

    MPI_Send(&N,1,MPI_INT,0,100,MPI_COMM_WORLD);
    MPI_Send(data,5*N+1,MPI_DOUBLE,0,100,MPI_COMM_WORLD);

    delete data;
}
#endif
