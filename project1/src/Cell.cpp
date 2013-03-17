#include "Cell.h"

#include <Atom.h>
#include <cstddef>
#include <inlines.h>
#include <iostream>
#include <System.h>
#include <threadcontrol.h>

Cell::Cell(System *system_)
{
    awesome_number = 1337;
    system = system_;
    num_atoms = 0;
    num_atoms_stored = 0;
    first_atom = NULL;
    reset();
 }

double dr_2, dr_6, dr_12, f, potential_energy, dr_12_inv, dr_6_inv;
double dx,dy,dz;

inline void Cell::calculate_force_between_atoms(Atom *atom0, Atom *atom1) {
    dx = atom0->r[0] - atom1->r[0];// + displacement_vector[0];
    dy = atom0->r[1] - atom1->r[1];// + displacement_vector[1];
    dz = atom0->r[2] - atom1->r[2];// + displacement_vector[2];

    if(dx>system->Lx/2)  dx -= system->Lx;
    if(dx<-system->Lx/2) dx += system->Lx;

    if(dy>system->Ly/2)  dy -= system->Ly;
    if(dy<-system->Ly/2) dy += system->Ly;

    if(dz>system->Lz/2)  dz -= system->Lz;
    if(dz<-system->Lz/2) dz += system->Lz;

    dr_2 = dx*dx + dy*dy + dz*dz;

    dr_6 = pow(dr_2,3);
    dr_12 = pow(dr_6,2);
    dr_12_inv = 1.0/dr_12;
    dr_6_inv = 1.0/dr_6;

    f = 24*(2.0*dr_12_inv-dr_6_inv)/dr_2;

    // potential_energy = 4*(dr_12_inv - dr_6_inv);

    atom0->a[0] += dx*f;
    atom0->a[1] += dy*f;
    atom0->a[2] += dz*f;

    atom1->a[0] -= dx*f;
    atom1->a[1] -= dy*f;
    atom1->a[2] -= dz*f;
}

void Cell::calculate_forces(System *system) {
    Atom *atom0, *atom1;
    for(unsigned long c=0;c<cells.size();c++) {
        int my_neighbor_cell_id = cells[c];
        Cell *my_neighbor_cell = system->thread_control->all_cells[my_neighbor_cell_id];
        if(my_neighbor_cell->forces_are_calculated) continue;

        atom0 = my_neighbor_cell->first_atom;
        while(atom0 != NULL) {
            atom1 = first_atom; // This is in my cell
            while(atom1 != NULL) {
                calculate_force_between_atoms(atom0, atom1);
                atom1 = atom1->next;
            }

            atom0 = atom0->next;
        }
    }

    atom0 = first_atom;
    while(atom0 != NULL) {
        atom1 = atom0->next;
        while(atom1 != NULL) {
            calculate_force_between_atoms(atom0, atom1);
            atom1 = atom1->next;
        }
        atom0 = atom0->next;
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
                if(cell_index == this->index) continue;

                vec displacement = zeros<vec>(3,1);

                displacement[0] = system->Lx*( -(di_p < di+i) + (di_p > di+i) );
                displacement[1] = system->Ly*( -(dj_p < dj+j) + (dj_p > dj+j) );
                displacement[2] = system->Lz*( -(dk_p < dk+k) + (dk_p > dk+k) );
                if(find(cells.begin(), cells.end(), cell_index) == cells.end()) {
                    cells.push_back(cell_index);
                    displacement_vectors.push_back(displacement);
                }
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
