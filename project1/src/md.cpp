#include <iostream>
#include "math.h"
#include "time.h"
#include <armadillo>
#include <fstream>
#include "Atom.h"
#include "System.h"
#include "StatisticsSampler.h"
#include <thermostat.h>
#include <CIniFile.h>
#include <unitconverter.h>
#include <settings.h>

#include <mpi.h>

using namespace arma;
using namespace std;

int main(int args, char *argv[]) {
    int numprocs = 1, my_rank = 0;
    MPI_Init(&args,&argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    Settings *settings = new Settings("../md.ini");

    System *system = new System(my_rank, settings);

	for(int i=0;i<timesteps;i++) {
        system->step();
	}

    MPI_Finalize();

    return 0;
}
