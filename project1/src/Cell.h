#pragma once
#include <Atom.h>

class Cell
{
public:
    Atom *first_atom;
    Atom *last_atom;
    int number_of_atoms;

    Cell();
    void add_atom(Atom *atom);
    void remove_atom(Atom *atom);
    void reset();
};
