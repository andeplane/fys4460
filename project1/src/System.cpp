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

void calculate_force_between_atoms(Atom *atom0, Atom *atom1, double &P) {
    double dr_2, dr_6, dr_12, f, potential_energy;

    vec dr = atom0->distanceToAtom(atom1);
    dr_2 = dot(dr,dr);

    dr_6 = pow(dr_2,3);
    dr_12 = pow(dr_6,2);

    f = 24*(2.0/dr_12-1.0/dr_6)/dr_2;

    potential_energy = 4*(1.0/dr_12 - 1.0/dr_6);

    atom0->a += f*dr;
    atom0->potential_energy += potential_energy;
    atom1->a -= f*dr;

    P += f*norm(dr,2);
}

void System::calculateAccelerations() {
#ifdef MPI_ENABLED
    P = 0;
    if(rank==0) {
        // Reset all atom accelerations
        for(int n=0;n<N;n++) {
            atoms[n]->a.zeros();
            atoms[n]->potential_energy = 0;
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
        sort_cells();
    }

    // time_t t0 = clock();
    // cout << "Time spent on sorting: " << ((double)clock()-t0)/CLOCKS_PER_SEC << endl;
    if(rank == 0) {
        for(int n=0;n<N;n++) {
            atoms[n]->v += atoms[n]->a*dt;
            atoms[n]->addR(atoms[n]->v*dt);
        }
    }

    calculateAccelerations();

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

#ifdef DEBUG_ME
        if(cell_index < 0 || cell_index >= cells.size()) {
            cout << "Problem with cell index for particle " << n << endl;
            cout << "Particle is at " << a->r << endl;
            cout << "which gives (i,j,k)=" << i << " " << j << " " << k << endl;
        }
#endif

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

#ifdef MPI_ENABLED
void System::send_particles_to_slaves() {
    sort_cells();

    double *positions_and_velocities = new double[6*N];
    int *indices = new int[N];
    int particles = 0;

    Cell *cell;

    for(int node_id=1;node_id<nodes;node_id++) {
        ThreadNode &node = thread_control->nodes[node_id];
        for(set<int>::iterator it=node.connected_cells.begin(); it!= node.connected_cells.end();it++) {
            int cell_index = *it;
            cell = cells[cell_index];
            for(int n=0;n<cell->atoms.size();n++) {
                positions_and_velocities[particles*6 + 0] = cell->atoms[n]->r(0);
                positions_and_velocities[particles*6 + 1] = cell->atoms[n]->r(1);
                positions_and_velocities[particles*6 + 2] = cell->atoms[n]->r(2);

                positions_and_velocities[particles*6 + 3] = cell->atoms[n]->v(0);
                positions_and_velocities[particles*6 + 4] = cell->atoms[n]->v(1);
                positions_and_velocities[particles*6 + 5] = cell->atoms[n]->v(2);
                indices[particles] = cell->atoms[n]->index;

                particles++;
            }
        }

        MPI_Send(&particles,1,MPI_INT,node_id,100,MPI_COMM_WORLD);
        MPI_Send(indices,particles,MPI_INT,node_id,100,MPI_COMM_WORLD);
        MPI_Send(positions_and_velocities,6*particles,MPI_DOUBLE,node_id,100,MPI_COMM_WORLD);
        particles = 0;
    }

    delete positions_and_velocities;
    delete indices;
}

void System::receive_particles_from_master() {
    for(int n=0;n<N;n++) {
        delete atoms[n];
    }

    atoms.clear();

    MPI_Status status;

    int particles = 0;
    MPI_Recv(&particles,1,MPI_INT,0,100,MPI_COMM_WORLD,&status);

    double *positions_and_velocities = new double[6*particles];
    int *indices = new int[particles];
    MPI_Recv(indices,particles,MPI_INT,0,100,MPI_COMM_WORLD,&status);
    MPI_Recv(positions_and_velocities,6*particles,MPI_DOUBLE,0,100,MPI_COMM_WORLD,&status);

    for(int n=0;n<particles;n++) {
        Atom *atom = new Atom(this);
        atom->r(0) = positions_and_velocities[6*n + 0];
        atom->r(1) = positions_and_velocities[6*n + 1];
        atom->r(2) = positions_and_velocities[6*n + 2];

        atom->v(0) = positions_and_velocities[6*n + 3];
        atom->v(1) = positions_and_velocities[6*n + 4];
        atom->v(2) = positions_and_velocities[6*n + 5];
        atom->index = indices[n];

        atoms.push_back(atom);
    }

    N = atoms.size();
    sort_cells();

    delete positions_and_velocities;
    delete indices;
}

void System::receive_particles_back_from_slaves() {
    MPI_Status status;

    for(int i=1;i<nodes;i++) {
        int particles = 0;
        double dP;
        MPI_Recv(&particles,1,MPI_INT,i,100,MPI_COMM_WORLD,&status);
        MPI_Recv(&dP,1,MPI_DOUBLE,i,100,MPI_COMM_WORLD,&status);
        P += dP;

        int *indices = new int[particles];
        double *accelerations = new double[4*particles];

        MPI_Recv(indices,particles,MPI_INT,i,100,MPI_COMM_WORLD,&status);
        MPI_Recv(accelerations,4*particles,MPI_DOUBLE,i,100,MPI_COMM_WORLD,&status);
        for(int j=0;j<particles;j++) {
            atoms[indices[j]]->a(0) += accelerations[4*j+0];
            atoms[indices[j]]->a(1) += accelerations[4*j+1];
            atoms[indices[j]]->a(2) += accelerations[4*j+2];
            atoms[indices[j]]->potential_energy += accelerations[4*j+3];
        }

        delete accelerations;
        delete indices;
    }

    P += N/V*T;
}

void System::send_particles_back_to_master() {
    int *indices = new int[N];
    double *accelerations = new double[4*N];

    for(int n=0;n<N;n++) {
        indices[n] = atoms[n]->index;
        accelerations[4*n+0] = atoms[n]->a(0);
        accelerations[4*n+1] = atoms[n]->a(1);
        accelerations[4*n+2] = atoms[n]->a(2);
        accelerations[4*n+3] = atoms[n]->potential_energy;
    }

    MPI_Send(&N,1,MPI_INT,0,100,MPI_COMM_WORLD);
    MPI_Send(&P,1,MPI_DOUBLE,0,100,MPI_COMM_WORLD);
    MPI_Send(indices,N,MPI_INT,0,100,MPI_COMM_WORLD);
    MPI_Send(accelerations,4*N,MPI_DOUBLE,0,100,MPI_COMM_WORLD);

    delete accelerations;
    delete indices;
}
#endif
