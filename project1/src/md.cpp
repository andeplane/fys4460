// #define ARMA_NO_DEBUG

#include <iostream>
#include "math.h"
#include "time.h"
#include <armadillo>
#include <fstream>
#include "Atom.h"
#include "System.h"
#include "StatisticsSampler.h"
#include <thermostat.h>
#include <mpi.h>

using namespace arma;
using namespace std;

int main(int args, char *argv[]) {


    int numprocs , my_rank;

    MPI_Init(&args,&argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int number_of_FCC_cells = 24;
    int T = 1;
    double dt = 0.01;
    int timesteps = 1000;

	// 108, 256, 500, 864, 1372
	// 2048, 2916, 4000, 5324, 6912
	// 8788, 10976, 13500
    System *system = new System(my_rank,numprocs,number_of_FCC_cells, T);

    MPI_Finalize();
    return 0;

    StatisticsSampler *sampler = NULL;
    ofstream *file = new ofstream;

    if(my_rank == 0) {
        sampler = new StatisticsSampler(system);

        file->open("pos.xyz");
        system->file = file;
        system->printPositionsToFile(file);
    }

    double t = 0;

    Thermostat thermostat(1.0,dt*20, dt);

	for(int i=0;i<timesteps;i++) {
        if(my_rank == 0 && !(i%(timesteps/100))) {
			printf("%d%%..",(100*i)/timesteps);
			fflush(stdout);
		}

		t+=dt;
        system->step(dt);

        if(my_rank == 0) {
            if(timesteps < 500) thermostat.apply(system->atoms);
            sampler->sample(t);

            system->printPositionsToFile(file);
        }

	}
	
    if(my_rank == 0) file->close();


    MPI_Finalize();

	return 0;
}
