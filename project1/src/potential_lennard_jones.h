#include <system.h>
#include <mdtimer.h>
#include <settings.h>

void System::calculate_accelerations() {
    mdtimer->start_forces();

    double dr2;
    double rr_cut = r_cut*r_cut;

    /* Reset the potential & forces */
    potential_energy = 0.0;
    pressure_forces = 0;
    memset(accelerations,0,num_atoms_local*3*sizeof(double));
    for (c=0; c<num_cells_including_ghosts_xyz; c++) head[c] = EMPTY;

    for (i=0; i<num_atoms_local+num_atoms_ghost; i++) {
        for (a=0; a<3; a++) mc[a] = (positions[i][a]+cell_length[a])/cell_length[a];

        cell_index_from_vector(mc,cell_index);

        // Set this atom at the head of the linked list
        linked_list[i] = head[cell_index];
        head[cell_index] = i;
    }

    double dr2_inverse, dr6_inverse;
    bool will_sample = settings->statistics_interval && steps % settings->statistics_interval == 0;
    double mass_inverse_24 = mass_inverse*24;

    for (mc[0]=1; mc[0]<=num_cells_local[0]; mc[0]++) {
        for (mc[1]=1; mc[1]<=num_cells_local[1]; mc[1]++) {
            for (mc[2]=1; mc[2]<=num_cells_local[2]; mc[2]++) {
                cell_index = mc[0]*num_cells_including_ghosts_yz+mc[1]*num_cells_including_ghosts[2]+mc[2];
                if ( head[cell_index] == EMPTY ) continue;

                // Loop through all neighbors of this cell. Note that i only sums over local atoms.
                for (mc1[0]=mc[0]-1; mc1[0]<=mc[0]+1; mc1[0]++) {
                    for (mc1[1]=mc[1]-1; mc1[1]<=mc[1]+1; mc1[1]++) {
                        for (mc1[2]=mc[2]-1; mc1[2]<=mc[2]+1; mc1[2]++) {
                            cell_index_2 = mc1[0]*num_cells_including_ghosts_yz+mc1[1]*num_cells_including_ghosts[2]+mc1[2];

                            if(head[cell_index_2] == EMPTY) continue;
                            i = head[cell_index];

                            while (i != EMPTY) {
                                j = head[cell_index_2];
                                while (j != EMPTY) {
#ifdef MANY_FROZEN_ATOMS
                                    if(i < j && !(atom_type[i]==FROZEN && atom_type[j]==FROZEN)) {
#else
                                    if(i < j) {
#endif
                                        /* Pair vector dr = r[i] - r[j] */
                                        dr[0] = positions[i][0]-positions[j][0];
                                        dr[1] = positions[i][1]-positions[j][1];
                                        dr[2] = positions[i][2]-positions[j][2];
                                        dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

                                        if (dr2<rr_cut) {
                                            bool is_local_atom = j < num_atoms_local;
                                            dr2_inverse = 1.0/dr2;
                                            dr6_inverse = dr2_inverse*dr2_inverse*dr2_inverse;

                                            double force = (2*dr6_inverse-1)*dr6_inverse*dr2_inverse*mass_inverse_24;

                                            if(will_sample) {
                                                double potential_energy_tmp = 4*dr6_inverse*(dr6_inverse - 1);
                                                if(is_local_atom) {
                                                    potential_energy += potential_energy_tmp;
                                                    pressure_forces += force*dr2;
                                                } else {
                                                    pressure_forces += 0.5*force*dr2;
                                                    potential_energy += 0.5*potential_energy_tmp;
                                                }
                                            }

                                            accelerations[3*i+0] += force*dr[0];
                                            accelerations[3*i+1] += force*dr[1];
                                            accelerations[3*i+2] += force*dr[2];

                                            if(is_local_atom) {
                                                accelerations[3*j+0] -= force*dr[0];
                                                accelerations[3*j+1] -= force*dr[1];
                                                accelerations[3*j+2] -= force*dr[2];
                                            }

                                        }
                                    } // if( i != j) {

                                    j = linked_list[j];
                                } // while (j != EMPTY) {
                                i = linked_list[i];

                            } // while (i != EMPTY) {


                            // cout << "Selected new i: " << i << endl;
                        } // for mc1[2]
                    } // for mc1[1]
                } // for mc1[0]
            } // for mc[2]
        } // for mc[1]
    } // for mc[0]
    // cout << "Atoms: " << potential_energy_count << endl;
    pressure_forces /= mass_inverse;
    mdtimer->end_forces();
}
