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
#define MAX_PARTICLE_NUM 1000000
#define MAX_CELL_NUM 1000
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
    int  *head;

    unsigned long steps;
    unsigned int myid;
    unsigned int node_index[3];
    unsigned int num_processors[3];
    unsigned int neighbor_nodes[6];
    short my_parity[3];
    double cell_length[3];
    double node_length[3];
    double system_length[3];
    int num_nodes;
    unsigned long num_atoms_local;
    unsigned long num_atoms_global;
    unsigned long num_atoms_ghost;
    long i,j,k,n,m,a,b,c, nx, ny, nz;

    double origo[3];
    double r_cut, dt, dt_half, potential_energy, t, volume;
    unsigned int mc[3];  // Usually cell index vector
    unsigned int mc1[3]; // Usually cell index vector
    unsigned int num_cells_including_ghosts_yz,cell_index, cell_index_2,num_cells_including_ghosts_xyz;
    unsigned int num_cells_local[3];
    unsigned int num_cells_including_ghosts[3];
    double dr[3];
    double shift_vector[6][3];
    double mpi_send_buffer[3*MAX_PARTICLE_NUM];
    double mpi_receive_buffer[3*MAX_PARTICLE_NUM];
    bool atom_moved[3*MAX_PARTICLE_NUM];
    double positions[3*MAX_PARTICLE_NUM];
    double accelerations[3*MAX_PARTICLE_NUM];
    double mass_inverse, pressure_forces;
    double velocities[3*MAX_PARTICLE_NUM];
    double initial_positions[3*MAX_PARTICLE_NUM];
    unsigned int move_queue[6][MAX_PARTICLE_NUM];
    bool frozen_atom[MAX_PARTICLE_NUM];
    int linked_list[MAX_PARTICLE_NUM];
    bool is_ghost_cell[MAX_CELL_NUM];


    System();
    void setup(int myid_, Settings *settings_);
    void step();
};
