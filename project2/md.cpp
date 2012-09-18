#include <iostream>
#include "math.h"
#include "time.h"
#include <armadillo>
#include <fstream>
#include "Atom.h"
#include "System.h"
#include "StatisticsSampler.h"

using namespace arma;
using namespace std;

int main(int argc, char *argv[]) {
	// 108, 256, 500, 864, 1372
	System *system = new System(108);

	StatisticsSampler *sampler = new StatisticsSampler(system);

	ofstream *file = new ofstream;
	file->open("pos.xyz");
	double dt = 0.005;
	double t = 0;

	system->printPositionsToFile(file);
	int timesteps = 1000;

	for(int i=0;i<timesteps;i++) {
		if(!(i%(timesteps/100))) {
			printf("%d%%..",(100*i)/timesteps);
			fflush(stdout);
		}

		t+=dt;
		system->step(dt);
		sampler->sample(t);

		system->printPositionsToFile(file);
	}
	
	file->close();

	return 0;
}