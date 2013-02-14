#include "Cell.h"
#include <cstddef>
#include <inlines.h>
#include <iostream>

Cell::Cell()
{
    initialized = 1337;

    reset();
    reset_atom_list();
 }


void Cell::calculate_force_between_atoms(Atom *atom0, Atom *atom1, double &P, const vec &displacement_vector) {
    double dr_2, dr_6, dr_12, f, potential_energy, dr_12_inv, dr_6_inv;
    double x,y,z;

    // atom0->distance_to_atom(atom1,displacement_vector,x,y,z);
    x = atom0->r(0) - atom1->r(0) + displacement_vector(0);
    y = atom0->r(1) - atom1->r(1) + displacement_vector(1);
    z = atom0->r(2) - atom1->r(2) + displacement_vector(2);

    dr_2 = x*x + y*y + z*z;

    dr_6 = pow(dr_2,3);
    dr_12 = pow(dr_6,2);
    dr_12_inv = 1.0/dr_12;
    dr_6_inv = 1.0/dr_6;


    f = 24*(2.0*dr_12_inv-dr_6_inv)/dr_2;

    potential_energy = 4*(dr_12_inv - dr_6_inv);

    atom0->a(0) += x*f;
    atom0->a(1) += y*f;
    atom0->a(2) += z*f;
    atom0->potential_energy += potential_energy;
    atom1->a(0) -= x*f;
    atom1->a(1) -= y*f;
    atom1->a(2) -= z*f;
    P += dr_2*f;
}

void Cell::calculate_forces(System *system) {
    Cell *cell;
    Atom *atom0, *atom1;
    const vec zero_vector = zeros<vec>(3,1);

    for(int c=0;c<cells.size();c++) {
        cell = system->cells[cells[c]];

        if(cell->forces_are_calculated) continue;

        const vec& displacement_vector = displacement_vectors[c];

        for(int k=0;k<cell->atoms.size();k++) {
            atom1 = cell->atoms[k];

            for(int i=0;i<atoms.size();i++) {
                atom0 = atoms[i];
                calculate_force_between_atoms(atom0, atom1, system->P, displacement_vector);
            }
        }
    }

    for(int i=0;i<atoms.size();i++) {
        atom0 = atoms[i];

        for(int j=i+1;j<atoms.size();j++) {

            atom1 = atoms[j];
            calculate_force_between_atoms(atom0, atom1, system->P, zero_vector);
        }
    }

    forces_are_calculated = true;
}

void Cell::find_neighbours(const int &c_x, const int &c_y, const int &c_z, System *system) {
    int di_p, dj_p, dk_p;

    for(int di=-1;di<=1;di++) {
        for(int dj=-1;dj<=1;dj++) {
            for(int dk=-1;dk<=1;dk++) {
                if(di == 0 && dj == 0 && dk == 0) continue;
                di_p = (i+di+10*c_x)%c_x;
                dj_p = (j+dj+10*c_y)%c_y;
                dk_p = (k+dk+10*c_z)%c_z;

                int cell_index = calculate_cell_index(di_p,dj_p,dk_p,c_x,c_y,c_z);

                vec displacement = zeros<vec>(3,1);

                displacement(0) = system->L*( -(di_p < di+i) + (di_p > di+i) );
                displacement(1) = system->L*( -(dj_p < dj+j) + (dj_p > dj+j) );
                displacement(2) = system->L*( -(dk_p < dk+k) + (dk_p > dk+k) );

                cells.push_back(cell_index);
                displacement_vectors.push_back(displacement);
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
