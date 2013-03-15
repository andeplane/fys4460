#pragma once

#include <vector>
#include <armadillo>
class System;
class Cell;
class Atom;

using namespace arma;
using namespace std;

class Cell
{
public:
    System *system;
    bool is_ghost_cell;
    bool forces_are_calculated;
    vector<int> cells;
    vector<vec> displacement_vectors;
    vector<Atom*> atoms;
    int i,j,k,index, node_id;

    Cell();
    void add_atom(Atom *atom);
    void remove_atom(Atom *atom);
    void calculate_force_between_atoms(Atom *atom0, Atom *atom1, double &P, const vec &displacement_vector);
    void reset();
    void reset_atom_list();
    void calculate_forces(System *system);
    void find_neighbours(const int &c_x, const int &c_y, const int &c_z, System *system);
};
