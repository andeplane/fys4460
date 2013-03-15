#pragma once
#include <vector>
#include <set>
#include <armadillo>

class System;
class Settings;
class Atom;
using namespace std;
using namespace arma;

class ThreadControl
{
public:
    int myid;
    int idx, idy, idz; // Geometric node index
    int num_nodes;     // Total number of nodes
    int num_atoms;
    int cells_x, cells_y, cells_z;
    vec3 origo;
    System *system;
    Settings *settings;

    vector<Cell*> my_cells;
    vector<Cell*> all_cells;
    vector<Cell*> ghost_cells;

    vector<Atom*> free_atoms;
    vector<Atom*> all_atoms;

    double *mpi_data;
    double *positions;
    double *velocities;
    double *initial_positions;
    int allocated_particle_data;

    ThreadControl();
    void setup(System *system);
    void setup_molecules();
    void setup_cells();
    void update_cells();
    void update_mpi();
    void update_local_cells();
    int cell_index_from_atom(Atom *atom);
    inline int cell_index_from_ijk(const int &i, const int &j, const int &k);
};
