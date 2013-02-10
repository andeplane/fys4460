#pragma once

#include <Atom.h>
#include <vector>
#include <armadillo>
#include <System.h>
class System;

using namespace std;
using namespace arma;

class Cell
{
public:
    Atom *first_atom;
    Atom *last_atom;
    int number_of_atoms;
    vector<int> cells;
    vector<Atom*> atoms;
    int i,j,k,index;
    bool forces_are_calculated;
    vec dr;
    int initialized;

    Cell();
    void add_atom(Atom *atom);
    void remove_atom(Atom *atom);
    void calculate_force_between_atoms(Atom *atom0, Atom *atom1, double &P);
    void reset();
    void reset_atom_list();
    void calculate_forces(System *system);
    void find_neighbours(const int c_x, const int c_y, const int c_z, System *system);
};
