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
    int numprocs = 1, myid = 0;
    MPI_Init(&args,&argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    cout << "We have " << numprocs << " processors, I am rank " << myid << endl;

    Settings *settings = new Settings("../md.ini");
    cout << "Settings loaded" << endl;
    System *system = new System(myid, settings);

    for(int i=0;i<settings->timesteps;i++) {
        system->step();
	}

    MPI_Finalize();

    return 0;
}
