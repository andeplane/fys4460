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
    vector<Cell*> cells;
    int i,j,k;
    bool forces_is_calculated;
    vec dr;

    Cell();
    void add_atom(Atom *atom);
    void remove_atom(Atom *atom);
    double calculate_force_between_atoms(Atom *atom0, Atom *atom1);
    void reset();
    void reset_atom_list();
    double calculate_forces();
    void find_neighbours(const int c_x, const int c_y, const int c_z, System *system);
};
