#include "Cell.h"
#include <cstddef>

Cell::Cell()
{
    reset();
}

void Cell::add_atom(Atom *atom) {
    if(first_atom == NULL) {
        first_atom = atom;
        last_atom = atom;

        atom->prev = NULL;
        atom->next = NULL;
    }
    else {
        last_atom->next = atom;
        last_atom = atom;

        atom->prev = last_atom;
        atom->next = NULL;
    }

    number_of_atoms++;
}

void Cell::remove_atom(Atom *atom) {
    if(first_atom == atom) {
        first_atom = atom->next;
    }
    else if(last_atom == atom) {
        last_atom = atom->prev;
    }

    if(atom->prev != NULL) atom->prev->next = atom->next;
    if(atom->next != NULL) atom->next->prev = atom->prev;

    atom->next = NULL;
    atom->prev = NULL;

    number_of_atoms--;
}

void Cell::reset() {
    first_atom = last_atom = NULL;
    number_of_atoms = 0;
}
