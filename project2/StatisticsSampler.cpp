#include "StatisticsSampler.h"
#include <armadillo>

using namespace arma;

StatisticsSampler::StatisticsSampler(System *system) {
	this->system = system;
	this->temperature = true;
	this->pressure = true;
	this->energy = true;

	this->temperatureFile = new ofstream("temperature.dat");
	this->pressureFile    = new ofstream("pressure.dat");
	this->energyFile      = new ofstream("energy.dat");
}

void StatisticsSampler::sample(double t) {
	this->calculateTemperature(t);
	this->calculateEnergy(t);
	this->calculatePressure(t);
}

void StatisticsSampler::calculateTemperature(double t) {
	if(!this->temperature) return;
	int N = this->system->getNumberOfAtoms();
	double vsquared = 0;

	Atom **atoms = this->system->getAtoms();
	double vmax = 0;
	for(int n=0;n<N;n++) {
		if(norm(atoms[n]->v,1) > vmax) vmax = norm(atoms[n]->v,1);

		vsquared += norm(atoms[n]->v,1);
	}

	vsquared/=(3*(N-1));

	*this->temperatureFile << t << " " << vsquared << endl;
}

void StatisticsSampler::calculateEnergy(double t) {
	if(!this->energy) return;
	int N = this->system->getNumberOfAtoms();
	double E = 0;

	Atom **atoms = this->system->getAtoms();
	double vmax = 0;
	for(int n=0;n<N;n++) {
		E += 0.5*atoms[n]->mass*norm(atoms[n]->v,1);
	}

	*this->energyFile << t << " " << E << endl;
}

void StatisticsSampler::calculatePressure(double t) {
	if(!this->pressure) return;

	
}

