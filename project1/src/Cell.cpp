#include "Cell.h"
#include <cstddef>
#include <inlines.h>

Cell::Cell()
{
    dr = zeros<vec>(3,1);
    reset_atom_list();
}

double Cell::calculate_force_between_atoms(Atom *atom0, Atom *atom1) {
    double dr_2, dr_6, dr_12, f, potential_energy;

    dr = atom0->distanceToAtom(atom1);

    dr_2 = dot(dr,dr);

    dr_6 = pow(dr_2,3);
    dr_12 = pow(dr_6,2);

    f = 24*(2.0/dr_12-1.0/dr_6)/dr_2;

    potential_energy = 2*(1.0/dr_12 - 1.0/dr_6);

    atom0->a += f*dr;
    atom0->potential_energy += potential_energy;
    atom1->a -= f*dr;
    atom1->potential_energy += potential_energy;
    return f*norm(dr,2);
}

double Cell::calculate_forces() {
    Cell *c;
    Atom *atom0, *atom1;
    double dP = 0;

    for(int i=0;i<cells.size();i++) {
        c = cells[i];
        if(c->forces_is_calculated) continue;
        Atom *atom0 = c->first_atom;
        while(atom0 != NULL) {
            atom1 = first_atom;
            while(atom1 != NULL) {

                dP += calculate_force_between_atoms(atom0, atom1);
                atom1 = atom1->next;
            }
            atom0 = atom0->next;
        }
    }

    atom0 = first_atom;
    while(atom0 != NULL) {
        atom1 = atom0->next;
        while(atom1 != NULL) {
            calculate_force_between_atoms(atom0,atom1);

            atom1 = atom1->next;
        }
        atom0 = atom0->next;
    }

    forces_is_calculated = true;
}

void Cell::find_neighbours(const int c_x, const int c_y, const int c_z, System *system) {

    Cell *c;
    for(int di=-1;di<=1;di++) {
        for(int dj=-1;dj<=1;dj++) {
            for(int dk=-1;dk<=1;dk++) {
                if(di == dj && di == dk && di == 0) continue;

                c = &system->cells[calculate_cell_index((i+di)%c_x,(j+dj)%c_y,((k+dk)%c_z),c_x,c_y,c_z)];
                cells.push_back(c);
            }
        }
    }
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
    forces_is_calculated = false;
}

void Cell::reset_atom_list() {
    first_atom = last_atom = NULL;
    number_of_atoms = 0;
}
