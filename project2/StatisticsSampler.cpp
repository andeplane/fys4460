#include "StatisticsSampler.h"
#include <armadillo>

using namespace arma;

StatisticsSampler::StatisticsSampler(System *system) {
	this->system = system;
	this->temperature = true;
	this->pressure = true;

	this->temperatureFile = new ofstream("temperature.dat");
	this->pressureFile    = new ofstream("pressure.dat");
}

void StatisticsSampler::sample(double t) {
	this->calculateTemperature(t);
}

void StatisticsSampler::calculateTemperature(double t) {
	if(!this->temperature) return;
	int N = this->system->getNumberOfAtoms();
	double vsquared = 0;

	Atom **atoms = this->system->getAtoms();

	for(int n=0;n<N;n++) {
		vsquared += norm(atoms[n]->v,1);
	}

	vsquared/=(3*(N-1));

	*this->temperatureFile << t << " " << vsquared << endl;
}

void StatisticsSampler::calculatePressure(double t) {
	if(!this->pressure) return;

	
}

