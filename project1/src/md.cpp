#include <iostream>
#include <math.h>
#include <time.h>
#include <armadillo>
#include <fstream>
#include <system.h>
#include <statisticssampler.h>
#include <thermostat.h>
#include <cinifile.h>
#include <unitconverter.h>
#include <settings.h>
#include <mpi.h>
#include <mdio.h>
#include <mdtimer.h>

using namespace arma;
using namespace std;

int main(int args, char *argv[]) {
    int numprocs = 1, myid = 0;
    MPI_Init(&args,&argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    Settings *settings = new Settings("../md.ini");
    System *system = new System();
    system->setup(myid, settings);

    for(int i=0;i<settings->timesteps;i++) {
        system->step();
        system->mdio->save_state_to_movie_file();
	}

    system->mdio->save_state_to_file_binary();

    if(myid==0) {
        cout << "System initialize : " << system->mdtimer->system_initialize << " s ( " << 100*system->mdtimer->fraction_system_initialize() << "%)" <<  endl;
        cout << "Disk IO           : " << system->mdtimer->io << " s ( " << 100*system->mdtimer->fraction_io() << "%)" <<  endl;
        cout << "Force calculation : " << system->mdtimer->forces << " s ( " << 100*system->mdtimer->fraction_forces() << "%)" <<  endl;
        cout << "Moving            : " << system->mdtimer->moving << " s ( " << 100*system->mdtimer->fraction_moving() << "%)" <<  endl;
        cout << "MPI communication : " << system->mdtimer->mpi << " s ( " << 100*system->mdtimer->fraction_mpi() << "%)" <<  endl;
    }

    MPI_Finalize();

    return 0;
}
