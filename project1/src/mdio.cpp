#include <mdio.h>
#include <system.h>
#include <settings.h>

MDIO::MDIO()
{

}

void MDIO::setup(System *system_) {
    system = system_;
    settings = system->settings;
    movie_file_open = false;
}

void MDIO::save_state_to_movie_file() {
    if(settings->create_movie && !(system->steps % settings->movie_every_n_frame)) {
        if(!movie_file_open) {
            char *filename = new char[100];
            sprintf(filename,"state_files/movie%04d.bin",system->myid);
            movie_file = new ofstream(filename,ios::out | ios::binary);
            movie_file_open = true;
            data = new double[3*settings->num_atoms_max];

            delete filename;
        }

        for(unsigned long n=0;n<system->num_atoms_local;n++) {
            data[3*n+0] = system->positions[3*n+0] + system->origo[0];
            data[3*n+1] = system->positions[3*n+1] + system->origo[1];
            data[3*n+2] = system->positions[3*n+2] + system->origo[2];
        }

        movie_file->write (reinterpret_cast<char*>(&system->num_atoms_local), sizeof(unsigned long));
        movie_file->write (reinterpret_cast<char*>(data), 3*system->num_atoms_local*sizeof(double));
        movie_file->flush();
    }
}

void MDIO::save_state_to_file_binary() {
//    if(system->myid==0) cout << "Saving state to file..." << endl;
//    ThreadControl *thread_control = system->thread_control;

//    char *filename = new char[100];
//    sprintf(filename,"state_files/state%04d.bin",system->myid);

//    ofstream file (filename, ios::out | ios::binary);
//    double *tmp_data = new double[9*thread_control->num_atoms];

//    int count = 0;
//    int atoms = 0;

//    for(unsigned int i=0;i<thread_control->my_cells.size();i++) {
//        Cell *cell = thread_control->my_cells[i];
//        Atom *atom = cell->first_atom;

//        while(atom != NULL) {
//            tmp_data[count++] = atom->r[0];
//            tmp_data[count++] = atom->r[1];
//            tmp_data[count++] = atom->r[2];

//            tmp_data[count++] = atom->r_initial[0];
//            tmp_data[count++] = atom->r_initial[1];
//            tmp_data[count++] = atom->r_initial[2];

//            tmp_data[count++] = atom->v[0];
//            tmp_data[count++] = atom->v[1];
//            tmp_data[count++] = atom->v[2];
//            atom = atom->next;
//            atoms++;
//        }
//    }

//    file.write (reinterpret_cast<char*>(&atoms), sizeof(int));
//    file.write (reinterpret_cast<char*>(tmp_data), 9*atoms*sizeof(double));

//    file.close();
//    delete tmp_data;
//    delete filename;
}

void MDIO::finalize() {
    if(movie_file_open) {
        movie_file->close();
    }

    if(system->myid != 0) return;
}
