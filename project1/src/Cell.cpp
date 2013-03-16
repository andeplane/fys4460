#include "Cell.h"

#include <Atom.h>
#include <cstddef>
#include <inlines.h>
#include <iostream>
#include <System.h>

Cell::Cell(System *system_)
{
    awesome_number = 1337;
    system = system_;
    num_atoms = 0;
    num_atoms_stored = 0;
    first_atom = NULL;
    reset();
 }

void Cell::calculate_force_between_atoms(Atom *atom0, Atom *atom1, const vec &displacement_vector) {
    double dr_2, dr_6, dr_12, f, potential_energy, dr_12_inv, dr_6_inv;
    double x,y,z;

    x = atom0->r[0] - atom1->r[0] + displacement_vector[0];
    y = atom0->r[1] - atom1->r[1] + displacement_vector[1];
    z = atom0->r[2] - atom1->r[2] + displacement_vector[2];

    dr_2 = x*x + y*y + z*z;

    dr_6 = pow(dr_2,3);
    dr_12 = pow(dr_6,2);
    dr_12_inv = 1.0/dr_12;
    dr_6_inv = 1.0/dr_6;

    f = 24*(2.0*dr_12_inv-dr_6_inv)/dr_2;

    potential_energy = 4*(dr_12_inv - dr_6_inv);

    atom0->a[0] += x*f;
    atom0->a[1] += y*f;
    atom0->a[2] += z*f;

    atom1->a[0] -= x*f;
    atom1->a[1] -= y*f;
    atom1->a[2] -= z*f;
}

void Cell::calculate_forces(System *system) {
    /*
    for(unsigned long c=0;c<cells.size();c++) {
        Cell *cell = system->cells[cells[c]];
        if(cell->forces_are_calculated) continue;

        const vec& displacement_vector = displacement_vectors[c];

        for(int k=0;k<cell->atoms.size();k++) {
            Atom *atom1 = cell->atoms[k];

            for(int i=0;i<atoms.size();i++) {
                Atom *atom0 = atoms[i];
                calculate_force_between_atoms(atom0, atom1, displacement_vector);
            }
        }
    }

    const vec zero_vector = zeros<vec>(3,1);

    for(int i=0;i<atoms.size();i++) {
        atom0 = atoms[i];

        for(int j=i+1;j<atoms.size();j++) {
            atom1 = atoms[j];
            calculate_force_between_atoms(atom0, atom1, zero_vector);
        }
    }

    forces_are_calculated = true;
    */
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

                displacement[0] = system->Lx*( -(di_p < di+i) + (di_p > di+i) );
                displacement[1] = system->Ly*( -(dj_p < dj+j) + (dj_p > dj+j) );
                displacement[2] = system->Lz*( -(dk_p < dk+k) + (dk_p > dk+k) );

                cells.push_back(cell_index);
                displacement_vectors.push_back(displacement);
            }
        }
    }

}

void Cell::add_atom(Atom *atom) {
    if(num_atoms_stored == 0) {
        first_atom = atom;
        last_atom = atom;
        atom->prev = NULL;
        atom->next = NULL;
        num_atoms_stored++;
        num_atoms++;
        atom->cell_index = this->index;
        return;
    }

    last_atom->next = atom;
    atom->prev = last_atom;
    last_atom = atom;
    atom->cell_index = this->index;

    num_atoms_stored++;
    num_atoms++;
}

void Cell::remove_atom(Atom *atom) {
    if(atom == first_atom && atom == last_atom) {
        first_atom = NULL;
        last_atom = NULL;
        atom->next = NULL;
        atom->prev = NULL;
        num_atoms--;
        num_atoms_stored--;
        return;
    }

    if(atom == first_atom) {
        if(atom->next) atom->next->prev = NULL;
        first_atom = atom->next;
        atom->next = NULL;
        atom->prev = NULL;
        num_atoms--;
        num_atoms_stored--;
        return;
    }
    if(atom == last_atom) {
        if(atom->prev) atom->prev->next = NULL;
        last_atom = atom->prev;
        atom->next = NULL;
        atom->prev = NULL;
        num_atoms--;
        num_atoms_stored--;
        return;
    }

    atom->prev->next = atom->next;
    if(atom->next) atom->next->prev = atom->prev;
    atom->next = NULL;
    atom->prev = NULL;
    num_atoms_stored--;
    num_atoms--;
}

void Cell::reset() {
    forces_are_calculated = false;
}
