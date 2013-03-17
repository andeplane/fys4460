#include <iostream>
#include "math.h"
#include "time.h"
#include <fstream>
#include "System.h"
#include <inlines.h>
#include <settings.h>
#include <mpi.h>
#include <mdio.h>
#include <unitconverter.h>

using namespace std;

System::System() {

}

void System::create_FCC() {
    UnitConverter *unit_converter = new UnitConverter();
    double xCell[4] = {0, 0.5, 0.5, 0};
    double yCell[4] = {0, 0.5, 0, 0.5};
    double zCell[4] = {0, 0, 0.5, 0.5};

    double rx,ry,rz;
    double T = unit_converter->temperature_from_SI(settings->temperature);

    for(int x = 0; x < settings->nodes_x*settings->unit_cells_x; x++) {
        for(int y = 0; y < settings->nodes_y*settings->unit_cells_y; y++) {
            for(int z = 0; z < settings->nodes_z*settings->unit_cells_z; z++) {
                for(int k = 0; k < 4; k++) {
                    // Set positions and type
                    rx = (x+xCell[k]) * settings->FCC_b - origo[0];
                    ry = (y+yCell[k]) * settings->FCC_b - origo[1];
                    rz = (z+zCell[k]) * settings->FCC_b - origo[2];

                    if(rx >= 0 && rx < Lx && ry >= 0 && ry < Ly && rz >= 0 && rz < Lz) {
                        positions[3*num_particles_local+0] = rx;
                        positions[3*num_particles_local+1] = ry;
                        positions[3*num_particles_local+2] = rz;

                        velocities[3*num_particles_local+0] = rnd->nextGauss()*sqrt(T);
                        velocities[3*num_particles_local+1] = rnd->nextGauss()*sqrt(T);
                        velocities[3*num_particles_local+2] = rnd->nextGauss()*sqrt(T);

                        num_particles_local++;
                    }
                }
            }
        }
    }

    MPI_Allreduce(&num_particles_local,&num_particles_global,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
}

void System::setup(int myid_, Settings *settings_) {
    myid = myid_;
    settings = settings_;
    dt = settings->dt;
    r_cut = settings->r_cut;
    num_particles_local = 0;
    num_particles_global = 0;
    num_particles_ghost = 0;

    positions = new double[3*settings->num_particles_max];
    accelerations = new double[3*settings->num_particles_max];
    velocities = new double[3*settings->num_particles_max];

    num_nodes = settings->nodes_x*settings->nodes_y*settings->nodes_z;

    steps = 0;
    rnd = new Random(-(myid+1));

    // Size of this node
    Lx = settings->unit_cells_x*settings->FCC_b;
    Ly = settings->unit_cells_y*settings->FCC_b;
    Lz = settings->unit_cells_z*settings->FCC_b;
    Lx_full = Lx*settings->nodes_x;
    Ly_full = Ly*settings->nodes_y;
    Lz_full = Lz*settings->nodes_z;

    node_idx = myid/(settings->nodes_y*settings->nodes_z);
    node_idy = (myid/settings->nodes_z) % settings->nodes_y;
    node_idz = myid%settings->nodes_z;

    origo[0] = (float)node_idx * Lx;
    origo[1] = (float)node_idy * Ly;
    origo[2] = (float)node_idz * Lz;

    mdio = new MDIO();
    mdio->setup(this);

    create_FCC();

    if(myid==0) cout << "System size: " << settings->nodes_x*Lx << " " << settings->nodes_y*Ly << " " << settings->nodes_z*Lz << endl;
    if(myid==0) cout << "Atoms: " << num_particles_global << endl;
}

double comm_time = 0;

void System::atom_copy() {

}

void System::atom_move() {

}

void System::half_kick() {

}

void System::kick() {

}

void System::calculate_accelerations() {

}

void System::step() {

}
