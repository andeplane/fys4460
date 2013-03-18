#pragma once

class Atom;
class ThreadControl;
class Settings;
class MDIO;
class MDTimer;
class Random;
class UnitConverter;

#include <fstream>
#include <vector>

#define EMPTY -1

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
    void init_parameters();
    void create_FCC();
    inline bool atom_did_change_node(double* ri, int ku);
    inline bool atom_should_be_copied(double *ri, int ku);
    inline void cell_index_from_ijk(const int &i, const int &j, const int &k, unsigned int &cell_index);
    inline void cell_index_from_vector(unsigned int *mc, unsigned int &cell_index);
public:
    Settings *settings;
    MDIO *mdio;
    Random *rnd;
    MDTimer *mdtimer;
    UnitConverter *unit_converter;

    unsigned int num_cells_including_ghosts_yz;
    unsigned int num_cells_including_ghosts_xyz;
    unsigned long steps;
    unsigned int myid;
    unsigned int node_index[3];
    unsigned int num_processors[3];
    unsigned int neighbor_nodes[6];
    double shift_vector[6][3];
    unsigned int mc[3];  // Usually cell index vector
    unsigned int mc1[3]; // Usually cell index vector
    double dr[3];
    short my_parity[3];
    double cell_length[3];
    double node_length[3];
    double system_length[3];
    unsigned int num_cells_local[3];
    unsigned int num_cells_including_ghosts[3];
    double *mpi_send_buffer;
    double *mpi_receive_buffer;
    bool *atom_moved;
    unsigned int *move_queue[6];
    int  *head;
    bool *is_ghost_cell;
    bool *frozen_atom;
    int *linked_list;

    int num_nodes;
    double mass_inverse, pressure_forces;
    unsigned long num_atoms_local;
    unsigned long num_atoms_global;
    unsigned long num_atoms_ghost;
    long i,j,k,n,m,a,b,c, nx, ny, nz;
    unsigned int cell_index, cell_index_2;

    double r_cut, dt, dt_half, potential_energy, t, volume;

    double *positions;
    double *initial_positions;
    double *accelerations;
    double *velocities;
    double origo[3];

    System();
    void setup(int myid_, Settings *settings_);
    void step();
};
