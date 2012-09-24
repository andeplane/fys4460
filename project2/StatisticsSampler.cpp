#include "StatisticsSampler.h"
#include <armadillo>

using namespace arma;

static int steps = 0;

StatisticsSampler::StatisticsSampler(System *system) {
	this->system = system;
	this->temperature = true;
	this->pressure = true;
	this->energy = true;
	this->printVelocities = true;

	this->temperatureFile = new ofstream("temperature.dat");
	this->pressureFile    = new ofstream("pressure.dat");
	this->energyFile      = new ofstream("energy.dat");
	// this->velocityFile      = new ofstream("velocities.dat");
	this->velocityFile = fopen("velocities.dat","w");

}

void StatisticsSampler::sample(double t) {
	this->calculateTemperature(t);
	this->calculateEnergy(t);
	this->calculatePressure(t);
	this->calculateVelocities(t);

	steps++;
}

void StatisticsSampler::calculateTemperature(double t) {
	if(!this->temperature) return;
	int N = this->system->N;
	double vsquared = 0;

	Atom **atoms = this->system->atoms;
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
	int N = this->system->N;
	double E = 0;

	Atom **atoms = this->system->atoms;
	double vmax = 0;
	for(int n=0;n<N;n++) {
		E += 0.5*atoms[n]->mass*norm(atoms[n]->v,1);
	}

	*this->energyFile << t << " " << E << endl;
}

void StatisticsSampler::calculatePressure(double t) {
	if(!this->pressure) return;

	
}

void StatisticsSampler::calculateVelocities(double t) {
	if(!this->printVelocities || steps % 50) return;

	int N = this->system->N;
	
	Atom **atoms = this->system->atoms;
	
	if(!steps) fprintf(this->velocityFile, "%d %d %d\n",N,N,N);
	Atom *atom;
	double vsq = 0;

	for(int n=0;n<N;n++) {
		atom = atoms[n];
		// vsq = norm(atom->v,2);

		fprintf(this->velocityFile, "%f %f %f \n",atom->v(0),atom->v(1),atom->v(2));
	}

}

