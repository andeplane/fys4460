#include <threadcontrol.h>

#include <Cell.h>
#include <inlines.h>
#include <system.h>
#include <Atom.h>
#include <settings.h>
#include <unitconverter.h>

ThreadControl::ThreadControl()
{

}

void ThreadControl::setup(System *system) {
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
    mpi_data = new double[9*settings->max_particle_num];

    if(myid==0) cout << "Setting up cells..." << endl;
    setup_cells();
    if(myid==0) cout << "Setting up molecules..." << endl;
    setup_molecules();
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
    double rx,ry,rz;
    int atom_count = 0;

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
                        Atom *atom = new Atom(this);

                        vx = system->rnd->nextGauss()*sqrt(T);
                        vy = system->rnd->nextGauss()*sqrt(T);
                        vz = system->rnd->nextGauss()*sqrt(T);
                        atom->r = &positions[3*atom_count];
                        atom->a = &accelerations[3*atom_count];
                        atom->v = &velocities[3*atom_count];
                        atom->r_initial = &initial_positions[3*atom_count];

                        atom->set_position(rx,ry,rz);
                        atom->set_initial_position(rx,ry,rz);
                        atom->set_velocity( vx , vy, vz );
                        atom->type = k>0; // For visualization

                        cell->add_atom(atom);
                        all_atoms.push_back(atom);
                        atom_count++;
                    }
                }
            }
        }
    }

    // TODO: Remove center of mass momentum
}

void ThreadControl::setup_cells() {
    int num_cells_global = num_nodes*settings->unit_cells_x*settings->unit_cells_y*settings->unit_cells_z;
    nodes_new_atoms_list.resize(num_nodes);
    dummy_cells.reserve(num_cells_global);

    for(int i=0;i<cells_x;i++) {
        for(int j=0;j<cells_y;j++) {
            for(int k=0;k<cells_z;k++) {
                int node_idx = i/settings->unit_cells_x;
                int node_idy = j/settings->unit_cells_y;
                int node_idz = k/settings->unit_cells_z;
                int node_id = node_idx*settings->nodes_y*settings->nodes_z + node_idy*settings->nodes_z + node_idz;

                Cell *c = new Cell(system);
                c->node_id = node_id;
                c->index = cell_index_from_ijk(i,j,k);
                cells.push_back(c);
                if(node_id == myid) {
                    my_cells.push_back(c);
                    all_cells.push_back(c);
                    c->is_ghost_cell = false;

                } else if(abs(node_idx - this->idx) <= 1 && abs(node_idy - this->idy) <= 1 && abs(node_idz - this->idz) <= 1) {
                    ghost_cells.push_back(c);
                    all_cells.push_back(c);
                    c->is_ghost_cell = true;

                } else {
                    all_cells.push_back(c);
                    c->is_ghost_cell = false;
                }
            }
        }
    }

    for(int i=0;i<my_cells.size();i++) {
        my_cells[i]->find_neighbours(cells_x,cells_y,cells_z, this);
    }
}
