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

    system->alpha = alpha;
    system->steepest_decent_at = steepest_decent_at;

    StatisticsSampler *sampler = NULL;
    ofstream *file = new ofstream;

    if(my_rank == 0) {
        sampler = new StatisticsSampler(system);
    }

    double t = 0;

    // Thermostat thermostat(1.0,dt*20, dt);

	for(int i=0;i<timesteps;i++) {
		t+=dt;
        system->step(dt);

        if(my_rank == 0) {
            // if(timesteps < 500) thermostat.apply(system->atoms);
            sampler->sample(t);

            if(print_positions) {
                system->printPositionsToFile(file);
            }
        }
	}
	
    if(my_rank == 0) file->close();

    MPI_Finalize();

    return 0;
}
