#include <threadcontrol.h>

#include <Cell.h>
#include <inlines.h>
#include <system.h>
#include <Atom.h>
#include <settings.h>
#include <unitconverter.h>
#include <mpi.h>
#include <algorithm>

ThreadControl::ThreadControl()
{

}

void ThreadControl::setup(System *system_) {
    system = system_;
    settings = system->settings;
    num_nodes = settings->nodes_x*settings->nodes_y*settings->nodes_z;
    myid = system->myid;

    idx = myid/(settings->nodes_y*settings->nodes_z); // Node id in x-direction
    idy = (myid/settings->nodes_z) % settings->nodes_y; // Node id in y-direction
    idz = myid % settings->nodes_z; // Node id in z-direction

    cout << "I am node " << myid << " with vector indices: " << idx << " " << idy << " " << idz << endl;

    origo[0] = idx*(system->Lx/settings->nodes_x);
    origo[1] = idy*(system->Ly/settings->nodes_y);
    origo[2] = idz*(system->Lz/settings->nodes_z);

    cout << "I am node " << myid << " with origo: " << origo[0] << " " << origo[1] << " " << origo[2] << endl;

    positions = new double[3*settings->max_particle_num];
    accelerations = new double[3*settings->max_particle_num];
    velocities = new double[3*settings->max_particle_num];
    initial_positions = new double[3*settings->max_particle_num];

    mpi_data_receive = new double[9*settings->max_particle_num];
    mpi_data_send = new double[9*settings->max_particle_num];

    if(myid==0) cout << "Setting up cells..." << endl;
    setup_cells();
    if(myid==0) cout << "Setting up molecules..." << endl;
    setup_molecules();
    // if(myid==0) cout << "Updating ghost cells..." << endl;
    // update_ghost_cells();
}

inline int ThreadControl::cell_index_from_ijk(const int &i, const int &j, const int &k) {
    return i*cells_y*cells_z + j*cells_z + k;
}

int ThreadControl::cell_index_from_atom(Atom *atom) {
    int i = atom->r[0] / system->Lx * cells_x;
    int j = atom->r[1] / system->Ly * cells_y;
    int k = atom->r[2] / system->Lz * cells_z;
    return cell_index_from_ijk(i,j,k);
}

void ThreadControl::setup_molecules() {
    double b = 1.545;
    UnitConverter *unit_converter = new UnitConverter();

    double T = unit_converter->temperature_from_SI(settings->temperature);

    double xCell[4] = {0, 0.5, 0.5, 0};
    double yCell[4] = {0, 0.5, 0, 0.5};
    double zCell[4] = {0, 0, 0.5, 0.5};
    double rx,ry,rz,vx,vy,vz;
    num_atoms = 0;

    for(int x = 0; x < settings->unit_cells_x; x++) {
        for(int y = 0; y < settings->unit_cells_y; y++) {
            for(int z = 0; z < settings->unit_cells_z; z++) {
                for(int k = 0; k < 4; k++) {
                    rx = (x+xCell[k]) * b;
                    ry = (y+yCell[k]) * b;
                    rz = (z+zCell[k]) * b;

                    int cell_x = rx / system->Lx * cells_x;
                    int cell_y = ry / system->Ly * cells_y;
                    int cell_z = rz / system->Lz * cells_z;

                    Cell *cell = all_cells[cell_index_from_ijk(cell_x,cell_y,cell_z)];

                    if(cell->node_id == myid) {
                        // This is our atom
                        Atom *atom = new Atom(system);

                        vx = system->rnd->nextGauss()*sqrt(T);
                        vy = system->rnd->nextGauss()*sqrt(T);
                        vz = system->rnd->nextGauss()*sqrt(T);
                        atom->r = &positions[3*num_atoms];
                        atom->a = &accelerations[3*num_atoms];
                        atom->v = &velocities[3*num_atoms];
                        atom->r_initial = &initial_positions[3*num_atoms];

                        atom->set_position(rx,ry,rz);
                        atom->set_initial_position(rx,ry,rz);
                        atom->set_velocity( vx , vy, vz );

                        cell->add_atom(atom);
                        all_atoms.push_back(atom);
                        num_atoms++;
                    }
                }
            }
        }
    }

    // TODO: Remove center of mass momentum
}

void ThreadControl::setup_cells() {
    cells_x = settings->nodes_x*settings->unit_cells_x;
    cells_y = settings->nodes_y*settings->unit_cells_y;
    cells_z = settings->nodes_z*settings->unit_cells_z;
    node_ghost_cell_list.resize(num_nodes);

    for(int i=0;i<cells_x;i++) {
        for(int j=0;j<cells_y;j++) {
            for(int k=0;k<cells_z;k++) {
                int node_idx = i/settings->unit_cells_x;
                int node_idy = j/settings->unit_cells_y;
                int node_idz = k/settings->unit_cells_z;
                int node_id = node_idx*settings->nodes_y*settings->nodes_z + node_idy*settings->nodes_z + node_idz;

                Cell *c = new Cell(system);
                c->i = i;
                c->j = j;
                c->k = k;
                c->node_id = node_id;
                c->index = cell_index_from_ijk(i,j,k);
                all_cells.push_back(c);
                c->is_ghost_cell = false;

                if(node_id == myid) {
                    my_cells.push_back(c);

                } else if(abs(node_idx - this->idx) <= 1 && abs(node_idy - this->idy) <= 1 && abs(node_idz - this->idz) <= 1) {
                    ghost_cells.push_back(c);
                    c->is_ghost_cell = true;
                }
            }
        }
    }

    for(unsigned long i=0;i<all_cells.size();i++) {
        my_cells[i]->find_neighbours(cells_x,cells_y,cells_z, system);
    }

    // Create node neighbor list with corresponding cells of interest for each node
    for(int di=-1;di<=1;di++) {
        for(int dj=-1;dj<=1;dj++) {
            for(int dk=-1;dk<=1;dk++) {
                int node_idx = (idx + di + settings->nodes_x) % settings->nodes_x;
                int node_idy = (idy + dj + settings->nodes_y) % settings->nodes_y;
                int node_idz = (idz + dk + settings->nodes_z) % settings->nodes_z;
                int node_id = node_idx*settings->nodes_y*settings->nodes_z + node_idy*settings->nodes_z + node_idz;
                if(find(neighbor_nodes.begin(), neighbor_nodes.end(), node_id) == neighbor_nodes.end()) {
                    neighbor_nodes.push_back(node_id);
                }

                // Loop through all of my cells and see if one of their neighbor cells belong to a neighbor node
                for(unsigned long i=0;i<my_cells.size();i++) {
                    Cell *my_cell = my_cells[i];
                    for(unsigned long j=0;j<my_cell->cells.size();j++) {
                        Cell *neighbor_cell = all_cells[my_cell->cells[j]];
                        if(neighbor_cell->node_id == node_id) {
                            // He is our friend, we should give him information from the current cell
                            node_ghost_cell_list[node_id].push_back(my_cell);
                        }
                    }
                }
            }
        }
    }

    // std::sort (neighbor_nodes.begin(), neighbor_nodes.end());
}

void ThreadControl::update_ghost_cells() {
    MPI_Status status;
    for(unsigned long i=0;i<ghost_cells.size();i++) {
        Cell *ghost_cell = ghost_cells[i];
        // Make all these atoms free
        free_atoms.insert(free_atoms.end(),ghost_cell->atoms.begin(),ghost_cell->atoms.end());
        ghost_cell->atoms.clear();
        ghost_cell->num_atoms = 0;
    }

    for(unsigned long i=0;i<neighbor_nodes.size();i++) {
        int node_id = neighbor_nodes[i];
        int atoms_sent = 0;
        int atoms_received = 0;
        vector<Cell*> &cells = node_ghost_cell_list[node_id];

        for(unsigned long c=0;c<cells.size();c++) {
            Cell *cell = cells[c];
            for(int n=0;n<cell->num_atoms;n++) {
                Atom *atom = cell->atoms[n];
                mpi_data_send[9*atoms_sent + 0] = atom->r[0];
                mpi_data_send[9*atoms_sent + 1] = atom->r[1];
                mpi_data_send[9*atoms_sent + 2] = atom->r[2];
                mpi_data_send[9*atoms_sent + 3] = atom->v[0];
                mpi_data_send[9*atoms_sent + 4] = atom->v[1];
                mpi_data_send[9*atoms_sent + 5] = atom->v[2];
                mpi_data_send[9*atoms_sent + 6] = atom->r_initial[0];
                mpi_data_send[9*atoms_sent + 7] = atom->r_initial[1];
                mpi_data_send[9*atoms_sent + 8] = atom->r_initial[2];
                atoms_sent++;
            }
        }

        MPI_Sendrecv(mpi_data_send,atoms_sent*9,MPI_DOUBLE,node_id,100,mpi_data_receive,atoms_received,MPI_DOUBLE,myid,100,MPI_COMM_WORLD,&status);
        atoms_received /= 9; // Each atom has 9 doubles from mpi_sendrecv

        double rx,ry,rz,vx,vy,vz,rx_initial,ry_initial,rz_initial;
        int cell_index;

        for(int j=0;j<atoms_received;j++) {
            Atom *atom;
            if(free_atoms.size() > 0) {
                atom = free_atoms.back();
                free_atoms.pop_back();
            } else {
                atom = new Atom(system);
                atom->r = &positions[3*num_atoms];
                atom->a = &accelerations[3*num_atoms];
                atom->v = &velocities[3*num_atoms];
                atom->r_initial = &initial_positions[3*num_atoms];

                all_atoms.push_back(atom);
                num_atoms++;
            }

            rx = mpi_data_receive[9*j + 0];
            ry = mpi_data_receive[9*j + 1];
            rz = mpi_data_receive[9*j + 2];
            vx = mpi_data_receive[9*j + 3];
            vy = mpi_data_receive[9*j + 4];
            vz = mpi_data_receive[9*j + 5];
            rx_initial = mpi_data_receive[9*j + 6];
            ry_initial = mpi_data_receive[9*j + 7];
            rz_initial = mpi_data_receive[9*j + 8];

            atom->set_position(rx,ry,rz);
            atom->set_initial_position(rx_initial,ry_initial,rz_initial);
            atom->set_velocity(vx , vy, vz);
            cell_index = cell_index_from_atom(atom);
            ((Cell*)all_cells[cell_index])->add_atom(atom);
        }
    }
}
