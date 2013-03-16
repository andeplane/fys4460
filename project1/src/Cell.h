#pragma once

#include <vector>
#include <armadillo>
class System;
class Atom;

using namespace arma;
using namespace std;

class Cell
{
public:
    System *system;
    Atom *first_atom;
    Atom *last_atom;
    bool is_ghost_cell;
    bool forces_are_calculated;
    vector<int> cells;
    vector<vec> displacement_vectors;
    vector<Atom*> new_atoms; // To be distributed to other nodes
    int i,j,k,index, node_id, num_atoms_stored, num_atoms;

    Cell(System *system_);
    void add_atom(Atom *atom);
    void remove_atom(Atom *atom);
    void calculate_force_between_atoms(Atom *atom0, Atom *atom1, const vec &displacement_vector);
    void reset();
    void calculate_forces(System *system);
    void find_neighbours(const int &c_x, const int &c_y, const int &c_z, System *system);
};
