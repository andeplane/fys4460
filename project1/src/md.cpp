#include <iostream>
#include "math.h"
#include "time.h"
#include <armadillo>
#include <fstream>
#include "System.h"
#include "StatisticsSampler.h"
#include <thermostat.h>
#include <CIniFile.h>
#include <unitconverter.h>
#include <settings.h>
#include <mpi.h>
#include <mdio.h>

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

    MPI_Finalize();

    return 0;
}
