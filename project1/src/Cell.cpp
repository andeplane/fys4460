#include "Cell.h"
#include <cstddef>
#include <inlines.h>
#include <iostream>

Cell::Cell()
{
    initialized = 1337;
    dr = zeros<vec>(3,1);
    reset();
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
    // atom1->potential_energy += potential_energy;

    return f*norm(dr,2);
}

double Cell::calculate_forces(System *system) {
    Cell *cell;
    Atom *atom0, *atom1;
    double dP = 0;
    for(int i=0;i<atoms.size();i++) {
        atom0 = atoms[i];

        for(int c=0;c<cells.size();c++) {
            cell = system->cells[cells[c]];

            if(cell->forces_are_calculated) continue;

            for(int k=0;k<cell->atoms.size();k++) {
                atom1 = cell->atoms[k];

                dP += calculate_force_between_atoms(atom0, atom1);
            }
        }
    }

    // Calculate force between atoms in this cell
    for(int i=0;i<atoms.size();i++) {
        atom0 = atoms[i];
        for(int j=i+1;j<atoms.size();j++) {

            atom1 = atoms[j];
            dP += calculate_force_between_atoms(atom0, atom1);
        }
    }

    forces_are_calculated = true;

    return dP;
}

void Cell::find_neighbours(const int c_x, const int c_y, const int c_z, System *system) {
    Cell *c;
    for(int di=-1;di<=1;di++) {
        for(int dj=-1;dj<=1;dj++) {
            for(int dk=-1;dk<=1;dk++) {
                if(di == 0 && dj == 0 && dk == 0) continue;
                int cell_index = calculate_cell_index((i+di+10*c_x)%c_x,(j+dj+10*c_y)%c_y,((k+dk+10*c_z)%c_z),c_x,c_y,c_z);

                cells.push_back(cell_index);
            }
        }
    }
}

void Cell::add_atom(Atom *atom) {
    atoms.push_back(atom);
    number_of_atoms++;
    /*
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
    */
}

void Cell::remove_atom(Atom *atom) {
    /*
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
    */
}

void Cell::reset() {
    forces_are_calculated = false;
}

void Cell::reset_atom_list() {
    atoms.clear();
    number_of_atoms = 0;
    /*
    first_atom = last_atom = NULL;
    number_of_atoms = 0;
    */
}
