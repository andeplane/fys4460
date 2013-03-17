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

    origo[0] = idx*(system->Lx/settings->nodes_x);
    origo[1] = idy*(system->Ly/settings->nodes_y);
    origo[2] = idz*(system->Lz/settings->nodes_z);

    positions = new double[3*settings->max_particle_num];
    accelerations = new double[3*settings->max_particle_num];
    velocities = new double[3*settings->max_particle_num];
    initial_positions = new double[3*settings->max_particle_num];

    mpi_particles_receive = new double[9*settings->max_particle_num];
    mpi_particles_send = new double[9*settings->max_particle_num];
    mpi_cells_send = new int[100000];
    mpi_cells_receive = new int[100000];
    setup_cells();
    setup_molecules();
    update_ghost_cells();
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
    int total_unit_cells_x = settings->unit_cells_x*settings->nodes_x;
    int total_unit_cells_y = settings->unit_cells_y*settings->nodes_y;
    int total_unit_cells_z = settings->unit_cells_z*settings->nodes_z;

    for(int x = 0; x < total_unit_cells_x; x++) {
        for(int y = 0; y < total_unit_cells_y; y++) {
            for(int z = 0; z < total_unit_cells_z; z++) {
                for(int k = 0; k < 4; k++) {

                    rx = (x+xCell[k]) * b + 0.001;
                    ry = (y+yCell[k]) * b + 0.001;
                    rz = (z+zCell[k]) * b + 0.001;

                    int cell_x = rx / system->Lx * cells_x;
                    int cell_y = ry / system->Ly * cells_y;
                    int cell_z = rz / system->Lz * cells_z;

                    Cell *cell = all_cells[cell_index_from_ijk(cell_x,cell_y,cell_z)];

                    if(cell->node_id == myid) {
                        // This is our atom
                        Atom *atom = create_new_atom();

                        vx = system->rnd->nextGauss()*sqrt(T);
                        vy = system->rnd->nextGauss()*sqrt(T);
                        vz = system->rnd->nextGauss()*sqrt(T);

                        atom->set_position(rx,ry,rz);
                        atom->set_initial_position(rx,ry,rz);
                        atom->set_velocity( vx , vy, vz );

                        cell->add_atom(atom);
                    }
                }
            }
        }
    }

    // TODO: Remove center of mass momentum
}

void ThreadControl::setup_cells() {
    cells_x = system->Lx / 3;
    cells_y = system->Ly / 3;
    cells_z = system->Lz / 3;

    int num_cells = cells_x*cells_y*cells_z;
    if(num_cells == 0) {
        cout << "Warning, system size is too small. Aborting!" << endl;
        MPI_Finalize();
        exit(0);
    }

    node_ghost_cell_list.resize(num_nodes);
    node_cell_list.resize(num_nodes);

    for(int i=0;i<cells_x;i++) {
        for(int j=0;j<cells_y;j++) {
            for(int k=0;k<cells_z;k++) {
                int node_idx = (float)i/cells_x*settings->nodes_x;
                int node_idy = (float)j/cells_y*settings->nodes_y;
                int node_idz = (float)k/cells_z*settings->nodes_z;

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
                    node_cell_list[node_id].push_back(c);
                    ghost_cells.push_back(c);
                    c->is_ghost_cell = true;
                }
            }
        }
    }

    for(unsigned long i=0;i<all_cells.size();i++) {
        Cell *cell = all_cells[i];
        cell->find_neighbours(cells_x,cells_y,cells_z, system);
    }

    // Create node neighbor list with corresponding cells of interest for each node
    for(int di=-1;di<=1;di++) {
        for(int dj=-1;dj<=1;dj++) {
            for(int dk=-1;dk<=1;dk++) {
                int node_idx = (idx + di + settings->nodes_x) % settings->nodes_x;
                int node_idy = (idy + dj + settings->nodes_y) % settings->nodes_y;
                int node_idz = (idz + dk + settings->nodes_z) % settings->nodes_z;
                int node_id = node_idx*settings->nodes_y*settings->nodes_z + node_idy*settings->nodes_z + node_idz;
                if(node_id == myid) continue;

                if(find(neighbor_nodes.begin(), neighbor_nodes.end(), node_id) == neighbor_nodes.end()) {
                    neighbor_nodes.push_back(node_id);
                }

                vector<Cell*> &this_node_ghost_cell_list = node_ghost_cell_list[node_id];

                // Loop through all of my cells and see if one of their neighbor cells belong to a neighbor node
                for(unsigned long i=0;i<my_cells.size();i++) {
                    Cell *my_cell = my_cells[i];
                    for(unsigned long j=0;j<my_cell->cells.size();j++) {
                        Cell *neighbor_cell = all_cells[my_cell->cells[j]];
                        if(neighbor_cell->node_id == node_id) {
                            // He is our friend, we should give him information from the current cell
                            if(find(this_node_ghost_cell_list.begin(), this_node_ghost_cell_list.end(), my_cell) == this_node_ghost_cell_list.end()) {
                                this_node_ghost_cell_list.push_back(my_cell);
                            }
                        }
                    }
                }
            }
        }
    }

    std::sort (neighbor_nodes.begin(), neighbor_nodes.end());
}

void ThreadControl::update_cells_local() {
    for(unsigned long i=0;i<my_cells.size();i++) {
        Cell *cell = my_cells[i];
        int atoms_in_this_cell = cell->new_atoms.size();

        for(int n=0;n<atoms_in_this_cell;n++) {
            Atom *atom = cell->new_atoms[n];
            Cell *old_cell = all_cells[atom->cell_index];
            old_cell->remove_atom(atom);
            cell->add_atom(atom);
        }
        cell->new_atoms.clear();
    }
}

void ThreadControl::update_cells_mpi() {
    MPI_Status status;

    for(unsigned long i=0;i<neighbor_nodes.size();i++) {
        int node_id = neighbor_nodes[i];
        if(node_id == myid) continue;

        int atoms_sent = 0;
        int atoms_received = 0;
        int cells_sent = 0;
        int cells_received = 0;

        vector<Cell*> &cells = node_cell_list[node_id];

        for(unsigned long c=0;c<cells.size();c++) {
            Cell *cell = cells[c];
            int particles_in_this_cell = cell->new_atoms.size();
            for(int n=0;n<particles_in_this_cell;n++) {
                Atom *atom = cell->new_atoms[n];
                Cell *old_cell = all_cells[atom->cell_index];

                mpi_particles_send[9*atoms_sent + 0] = atom->r[0];
                mpi_particles_send[9*atoms_sent + 1] = atom->r[1];
                mpi_particles_send[9*atoms_sent + 2] = atom->r[2];
                mpi_particles_send[9*atoms_sent + 3] = atom->v[0];
                mpi_particles_send[9*atoms_sent + 4] = atom->v[1];
                mpi_particles_send[9*atoms_sent + 5] = atom->v[2];
                mpi_particles_send[9*atoms_sent + 6] = atom->r_initial[0];
                mpi_particles_send[9*atoms_sent + 7] = atom->r_initial[1];
                mpi_particles_send[9*atoms_sent + 8] = atom->r_initial[2];
                atoms_sent++;

                old_cell->remove_atom(atom);
                free_atoms.push_back(atom);
            }

            cell->new_atoms.clear();
            mpi_cells_send[2*cells_sent+0] = cell->index;
            mpi_cells_send[2*cells_sent+1] = particles_in_this_cell;
            cells_sent++;
        }

        cells_sent *= 2; // Two values per cell <cell_index, num_particles>
        atoms_sent *= 9; // 9 doubles each

        if(node_id > myid) {
            MPI_Send(&cells_sent,1,MPI_INT,node_id,100,MPI_COMM_WORLD);
            MPI_Send(mpi_cells_send,cells_sent,MPI_INT,node_id,100,MPI_COMM_WORLD);

            MPI_Send(&atoms_sent,1,MPI_INT,node_id,100,MPI_COMM_WORLD);
            MPI_Send(mpi_particles_send,atoms_sent,MPI_DOUBLE,node_id,100,MPI_COMM_WORLD);

            MPI_Recv(&cells_received,1,MPI_INT,node_id,100,MPI_COMM_WORLD,&status);
            MPI_Recv(mpi_cells_receive,cells_received,MPI_INT,node_id,100,MPI_COMM_WORLD,&status);

            MPI_Recv(&atoms_received,1,MPI_INT,node_id,100,MPI_COMM_WORLD,&status);
            MPI_Recv(mpi_particles_receive,atoms_received,MPI_DOUBLE,node_id,100,MPI_COMM_WORLD,&status);
        } else {
            MPI_Recv(&cells_received,1,MPI_INT,node_id,100,MPI_COMM_WORLD,&status);
            MPI_Recv(mpi_cells_receive,cells_received,MPI_INT,node_id,100,MPI_COMM_WORLD,&status);

            MPI_Recv(&atoms_received,1,MPI_INT,node_id,100,MPI_COMM_WORLD,&status);
            MPI_Recv(mpi_particles_receive,atoms_received,MPI_DOUBLE,node_id,100,MPI_COMM_WORLD,&status);

            MPI_Send(&cells_sent,1,MPI_INT,node_id,100,MPI_COMM_WORLD);
            MPI_Send(mpi_cells_send,cells_sent,MPI_INT,node_id,100,MPI_COMM_WORLD);

            MPI_Send(&atoms_sent,1,MPI_INT,node_id,100,MPI_COMM_WORLD);
            MPI_Send(mpi_particles_send,atoms_sent,MPI_DOUBLE,node_id,100,MPI_COMM_WORLD);
        }

        cells_received /= 2; // 2 ints per cell

        double rx,ry,rz,vx,vy,vz,rx_initial,ry_initial,rz_initial;
        int new_particles = 0;

        for(int j=0;j<cells_received;j++) {
            int cell_index = mpi_cells_receive[2*j+0];
            int new_particles_in_this_cell = mpi_cells_receive[2*j+1];

            Cell *cell = all_cells[cell_index];

            for(int k=0;k<new_particles_in_this_cell;k++) {
                Atom *atom;
                if(free_atoms.size() > 0) {
                    atom = free_atoms.back();
                    free_atoms.pop_back();
                } else {
                    atom = create_new_atom();
                }

                cell->add_atom(atom);

                rx 		   = mpi_particles_receive[9*new_particles + 0];
                ry 		   = mpi_particles_receive[9*new_particles + 1];
                rz 		   = mpi_particles_receive[9*new_particles + 2];
                vx 		   = mpi_particles_receive[9*new_particles + 3];
                vy 		   = mpi_particles_receive[9*new_particles + 4];
                vz 		   = mpi_particles_receive[9*new_particles + 5];
                rx_initial = mpi_particles_receive[9*new_particles + 6];
                ry_initial = mpi_particles_receive[9*new_particles + 7];
                rz_initial = mpi_particles_receive[9*new_particles + 8];

                atom->set_position(rx,ry,rz);
                atom->set_initial_position(rx_initial,ry_initial,rz_initial);
                atom->set_velocity(vx , vy, vz);
                atom->set_acceleration(0,0,0);

                new_particles++;
            }
        }
    }
}

void ThreadControl::update_ghost_cells() {
    MPI_Status status;

    for(unsigned long i=0;i<neighbor_nodes.size();i++) {
        int node_id = neighbor_nodes[i];
        if(node_id == myid) continue;

        int atoms_sent = 0;
        int atoms_received = 0;
        int cells_sent = 0;
        int cells_received = 0;

        vector<Cell*> &cells = node_ghost_cell_list[node_id];

        for(unsigned long c=0;c<cells.size();c++) {
            Cell *cell = cells[c];
            int particles_in_this_cell = 0;
            Atom *atom = cell->first_atom;
            while(atom != NULL) {
                mpi_particles_send[9*atoms_sent + 0] = atom->r[0];
                mpi_particles_send[9*atoms_sent + 1] = atom->r[1];
                mpi_particles_send[9*atoms_sent + 2] = atom->r[2];
                mpi_particles_send[9*atoms_sent + 3] = atom->v[0];
                mpi_particles_send[9*atoms_sent + 4] = atom->v[1];
                mpi_particles_send[9*atoms_sent + 5] = atom->v[2];
                mpi_particles_send[9*atoms_sent + 6] = atom->r_initial[0];
                mpi_particles_send[9*atoms_sent + 7] = atom->r_initial[1];
                mpi_particles_send[9*atoms_sent + 8] = atom->r_initial[2];
                particles_in_this_cell++;
                atoms_sent++;

                atom = atom->next;
            }

            mpi_cells_send[2*cells_sent+0] = cell->index;
            mpi_cells_send[2*cells_sent+1] = particles_in_this_cell;
            cells_sent++;
        }

        cells_sent *= 2; // Two values per cell <cell_index, num_particles>
        atoms_sent *= 9; // 9 doubles each

        if(node_id > myid) {
            MPI_Send(&cells_sent,1,MPI_INT,node_id,100,MPI_COMM_WORLD);
            MPI_Send(mpi_cells_send,cells_sent,MPI_INT,node_id,100,MPI_COMM_WORLD);

            MPI_Send(&atoms_sent,1,MPI_INT,node_id,100,MPI_COMM_WORLD);
            MPI_Send(mpi_particles_send,atoms_sent,MPI_DOUBLE,node_id,100,MPI_COMM_WORLD);

            MPI_Recv(&cells_received,1,MPI_INT,node_id,100,MPI_COMM_WORLD,&status);
            MPI_Recv(mpi_cells_receive,cells_received,MPI_INT,node_id,100,MPI_COMM_WORLD,&status);

            MPI_Recv(&atoms_received,1,MPI_INT,node_id,100,MPI_COMM_WORLD,&status);
            MPI_Recv(mpi_particles_receive,atoms_received,MPI_DOUBLE,node_id,100,MPI_COMM_WORLD,&status);
        } else {
            MPI_Recv(&cells_received,1,MPI_INT,node_id,100,MPI_COMM_WORLD,&status);
            MPI_Recv(mpi_cells_receive,cells_received,MPI_INT,node_id,100,MPI_COMM_WORLD,&status);

            MPI_Recv(&atoms_received,1,MPI_INT,node_id,100,MPI_COMM_WORLD,&status);
            MPI_Recv(mpi_particles_receive,atoms_received,MPI_DOUBLE,node_id,100,MPI_COMM_WORLD,&status);

            MPI_Send(&cells_sent,1,MPI_INT,node_id,100,MPI_COMM_WORLD);
            MPI_Send(mpi_cells_send,cells_sent,MPI_INT,node_id,100,MPI_COMM_WORLD);

            MPI_Send(&atoms_sent,1,MPI_INT,node_id,100,MPI_COMM_WORLD);
            MPI_Send(mpi_particles_send,atoms_sent,MPI_DOUBLE,node_id,100,MPI_COMM_WORLD);
        }

        cells_received /= 2; // Two values per cell <cell_index, num_particles>

        double rx,ry,rz,vx,vy,vz,rx_initial,ry_initial,rz_initial;
        int new_particles = 0;

        for(int j=0;j<cells_received;j++) {
            int cell_index = mpi_cells_receive[2*j+0];
            int particles_in_this_cell = mpi_cells_receive[2*j+1];

            Cell *cell = all_cells[cell_index];
            Atom *atom = cell->first_atom;

            for(int k=0;k<particles_in_this_cell;k++) {
                if(atom == NULL) {
                    // Cell doesn't have any free atoms
                    if(free_atoms.size() > 0) {
                        atom = free_atoms.back();
                        free_atoms.pop_back();
                    } else {
                        atom = create_new_atom();
                    }

                    cell->add_atom(atom);
                }

                rx 		   = mpi_particles_receive[9*new_particles + 0];
                ry 		   = mpi_particles_receive[9*new_particles + 1];
                rz 		   = mpi_particles_receive[9*new_particles + 2];
                vx 		   = mpi_particles_receive[9*new_particles + 3];
                vy 		   = mpi_particles_receive[9*new_particles + 4];
                vz 		   = mpi_particles_receive[9*new_particles + 5];
                rx_initial = mpi_particles_receive[9*new_particles + 6];
                ry_initial = mpi_particles_receive[9*new_particles + 7];
                rz_initial = mpi_particles_receive[9*new_particles + 8];

                atom->set_position(rx,ry,rz);
                atom->set_initial_position(rx_initial,ry_initial,rz_initial);
                atom->set_velocity(vx , vy, vz);
                atom->set_acceleration(0,0,0);

                new_particles++;
                atom = atom->next;
            }
        }
    }
}

Atom *ThreadControl::create_new_atom() {
    Atom *atom = new Atom(system);
    atom->r = &positions[3*num_atoms];
    atom->a = &accelerations[3*num_atoms];
    atom->v = &velocities[3*num_atoms];
    atom->r_initial = &initial_positions[3*num_atoms];

    atom->index = num_atoms;
    num_atoms++;

    atom->set_acceleration(0,0,0);
    all_atoms.push_back(atom);

    return atom;
}

void ThreadControl::reset_forces() {
    for(unsigned long i=0;i<all_cells.size();i++) {
        Cell *cell = all_cells[i];
        cell->forces_are_calculated = false;

        Atom *atom = cell->first_atom;
        while(atom != NULL) {
            atom->a[0] = 0;
            atom->a[1] = 0;
            atom->a[2] = 0;
            atom = atom->next;
        }
    }
}
