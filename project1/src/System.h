#pragma once

class Atom;
class ThreadControl;
class Settings;
class MDIO;

#include <fstream>
#include <random.h>
#include <vector>

using namespace std;

class System {
private:
    void initialize();
    void calculate_accelerations();
    void half_kick();
    void full_kick();
    void mpi_move();
    void mpi_copy();
    void move();
    void set_topology();
    void create_FCC();
    inline int atom_did_change_node(double* ri, int ku);
public:
    Settings *settings;
    MDIO *mdio;
    Random *rnd;

    long steps;
    int myid;
    int node_index[3];
    int num_processors[3];
    int neighbor_nodes[6];
    double shift_vector[6][3];
    short my_parity[3];
    double box_length[3];
    double box_length_full[3];
    double *mpi_send_buffer;
    double *mpi_receive_buffer;
    bool *atom_moved;
    int *move_queue[6];

    int num_nodes;
    unsigned long num_atoms_local;
    unsigned long num_atoms_global;
    unsigned long num_atoms_ghost;
    unsigned long i,j,k,n,m,a,b,c, nx, ny, nz;

    double r_cut, dt, dt_half;

    double *positions;
    double *accelerations;
    double *velocities;
    double origo[3];

    System();
    void setup(int myid_, Settings *settings_);
    void step();
};
