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

System::System(int rank_, int nodes_, int number_of_FCC_cells_, double T_) {
    number_of_FCC_cells = number_of_FCC_cells_;
    T = T_;
    rank = rank_;
    nodes = nodes_;
    initialize();
}

void System::calculateAccelerations() {
    P = N/V*T;
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

void System::send_particles_to_slaves() {
    sort_cells();
    int cell_rows_per_node = cells_z/(nodes-1);
    int current_node = 1;
    int cell_rows_so_far = 0;
    int cell_index, particles;
    int max_node_index = nodes - 1;
    double *positions_and_velocities = new double[6*N];
    int *indices = new int[N];

    printf("We have %d particles on %d nodes. We have %d cells in z-direction\n",N,nodes-1,cells_z);

    Cell *cell;
    particles = 0;

    for(int k=0;k<cells_z;k++) {
        cell_rows_so_far++;
        for(int i=0;i<cells_x;i++) {
            for(int j=0;j<cells_y;j++) {
                cell_index = calculate_cell_index(i,j,k,cells_x,cells_y,cells_z);
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
        }

        if(cell_rows_so_far >= cell_rows_per_node && current_node < max_node_index) {
            // Maximum cell rows for this node, send data
            MPI_Send(&particles,1,MPI_INT,current_node,100,MPI_COMM_WORLD);
            MPI_Send(indices,particles,MPI_INT,current_node,100,MPI_COMM_WORLD);
            MPI_Send(positions_and_velocities,6*particles,MPI_DOUBLE,current_node,100,MPI_COMM_WORLD);
            current_node++;
            cell_rows_so_far = 0;
            particles = 0;
        }
    }

    if(particles > 0) {
        // Send the rest
        MPI_Send(&particles,1,MPI_INT,current_node,100,MPI_COMM_WORLD);

        MPI_Send(indices,particles,MPI_INT,current_node,100,MPI_COMM_WORLD);
        MPI_Send(positions_and_velocities,6*particles,MPI_DOUBLE,current_node,100,MPI_COMM_WORLD);
    }
}

void System::receive_particles_from_master() {
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
    cout << "I am " << rank << " with " << N << "particles." << endl;
    sort_cells();
}

void System::receive_particles_back_from_slaves() {
    MPI_Status status;

    for(int i=1;i<nodes;i++) {
        int particles = 0;

        MPI_Recv(&particles,1,MPI_INT,i,100,MPI_COMM_WORLD,&status);

        int *indices = new int[particles];
        double *accelerations = new double[3*particles];

        MPI_Recv(indices,particles,MPI_INT,i,100,MPI_COMM_WORLD,&status);
        MPI_Recv(accelerations,3*particles,MPI_DOUBLE,i,100,MPI_COMM_WORLD,&status);
        for(int j=0;j<particles;j++) {
            atoms[indices[j]]->a(0) = accelerations[3*j+0];
            atoms[indices[j]]->a(1) = accelerations[3*j+1];
            atoms[indices[j]]->a(2) = accelerations[3*j+2];
        }

    }
}

void System::send_particles_back_to_master() {
    int *indices = new int[N];
    double *accelerations = new double[3*N];

    for(int n=0;n<N;n++) {
        indices[n] = atoms[n]->index;
        accelerations[3*n+0] = atoms[n]->a(0);
        accelerations[3*n+1] = atoms[n]->a(1);
        accelerations[3*n+2] = atoms[n]->a(2);
    }

    MPI_Send(&N,1,MPI_INT,0,100,MPI_COMM_WORLD);
    MPI_Send(indices,N,MPI_INT,0,100,MPI_COMM_WORLD);
    MPI_Send(accelerations,3*N,MPI_DOUBLE,0,100,MPI_COMM_WORLD);
}
