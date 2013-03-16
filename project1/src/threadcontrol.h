#pragma once
#include <vector>
#include <set>
#include <armadillo>

class System;
class Settings;
class Atom;
class Cell;
using namespace std;
using namespace arma;

class ThreadControl
{
public:
    System *system;
    Settings *settings;

    vector<Cell*> my_cells;
    vector<Cell*> all_cells;
    vector<Cell*> ghost_cells;
    vector< vector< Cell*> > node_ghost_cell_list;
    vector<int> neighbor_nodes;
    vector<Atom*> free_atoms;
    vector<Atom*> all_atoms;

    double *mpi_data_receive;
    double *mpi_data_send;
    double *positions;
    double *accelerations;
    double *velocities;
    double *initial_positions;

    int allocated_particle_data;
    int myid;
    int idx, idy, idz; // Geometric node index
    int num_nodes;     // Total number of nodes
    int num_atoms;
    int cells_x, cells_y, cells_z;
    vec3 origo;

    ThreadControl();
    void setup(System *system);
    void setup_molecules();
    void setup_cells();
    void update_cells();

    void update_ghost_cells();
    void update_local_cells();
    int cell_index_from_atom(Atom *atom);
    inline int cell_index_from_ijk(const int &i, const int &j, const int &k);
};
