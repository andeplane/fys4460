#include <mdio.h>
#include <System.h>
#include <settings.h>
#include <threadcontrol.h>
#include <Cell.h>
#include <Atom.h>

MDIO::MDIO()
{

}

void MDIO::setup(System *system_) {
    system = system_;
    settings = system->settings;
}

void MDIO::save_state_to_file_binary() {
    if(system->myid==0) cout << "Saving state to file..." << endl;
    ThreadControl *thread_control = system->thread_control;

    char *filename = new char[100];
    sprintf(filename,"state_files/state%04d.bin",system->myid);

    ofstream file (filename, ios::out | ios::binary);
    double *tmp_data = new double[9*thread_control->num_atoms];

    int count = 0;
    int atoms = 0;

    for(unsigned int i=0;i<thread_control->my_cells.size();i++) {
        Cell *cell = thread_control->my_cells[i];
        Atom *atom = cell->first_atom;

        while(atom != NULL) {
            tmp_data[count++] = atom->r[0];
            tmp_data[count++] = atom->r[1];
            tmp_data[count++] = atom->r[2];

            tmp_data[count++] = atom->r_initial[0];
            tmp_data[count++] = atom->r_initial[1];
            tmp_data[count++] = atom->r_initial[2];

            tmp_data[count++] = atom->v[0];
            tmp_data[count++] = atom->v[1];
            tmp_data[count++] = atom->v[2];
            atom = atom->next;
            atoms++;
        }
    }

    file.write (reinterpret_cast<char*>(&atoms), sizeof(int));
    file.write (reinterpret_cast<char*>(tmp_data), 9*atoms*sizeof(double));

    file.close();
    delete tmp_data;
    delete filename;
}
