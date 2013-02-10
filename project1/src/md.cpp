#include <iostream>
#include "math.h"
#include "time.h"
#include <armadillo>
#include <fstream>
#include "Atom.h"
#include "System.h"
#include "StatisticsSampler.h"
#include <thermostat.h>
#ifdef MPI_ENABLED
#include <mpi.h>
#endif

using namespace arma;
using namespace std;

int main(int args, char *argv[]) {
    // QTime myTimer;
    // myTimer.start();

    int numprocs = 1, my_rank = 0;
#ifdef MPI_ENABLED
    MPI_Init(&args,&argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

    int number_of_FCC_cells = 6;
    double T = 1.0;
    double dt = 0.01;
    int timesteps = 10;

    System *system = new System(my_rank,numprocs,dt,number_of_FCC_cells, T);

    StatisticsSampler *sampler = NULL;
    ofstream *file = new ofstream;

    if(my_rank == 0) {
        sampler = new StatisticsSampler(system);

        file->open("pos.xyz");
        system->file = file;
        // system->printPositionsToFile(file);
    }

    double t = 0;

    Thermostat thermostat(1.0,dt*20, dt);

	for(int i=0;i<timesteps;i++) {
		t+=dt;
        system->step(dt);

        if(my_rank == 0) {
            if(timesteps < 500) thermostat.apply(system->atoms);
            // sampler->sample(t);

            // system->printPositionsToFile(file);
        }

	}
	
    if(my_rank == 0) file->close();

#ifdef MPI_ENABLED
    MPI_Finalize();
#endif
    // do something..
    /*
    int nMilliseconds = myTimer.elapsed();
    cout << "Program used " << nMilliseconds/1000 << " seconds." << endl;
    */
	return 0;
}
