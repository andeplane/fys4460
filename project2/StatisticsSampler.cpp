#include "StatisticsSampler.h"
#include <armadillo>

using namespace arma;

static int steps = 0;

StatisticsSampler::StatisticsSampler(System *system) {
	this->system = system;
	this->temperature = true;
	this->pressure = true;
	this->energy = true;
	this->printVelocities = false;
	this->diffusionConstant = true;

	this->temperatureFile = fopen("temperature.dat","w");
	this->pressureFile    = fopen("pressure.dat","w");
	this->energyFile      = fopen("energy.dat","w");
	this->velocityFile = fopen("velocity.dat","w");
	this->diffusionFile = fopen("diffusion.dat","w");
}

void StatisticsSampler::sample(double t) {
	this->calculateTemperature(t);
	this->calculateEnergy(t);
	this->calculatePressure(t);
	this->calculateVelocities(t);
	this->calculateDiffusionConstant(t);

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

	fprintf(this->temperatureFile, "%f %f \n",t, vsquared);
}

void StatisticsSampler::calculateEnergy(double t) {
	if(!this->energy) return;
	int N = this->system->N;
	double E = 0, Ek=0,Ep=0, Ek_temp, Ep_temp;

	Atom **atoms = this->system->atoms;
	double vmax = 0;
	for(int n=0;n<N;n++) {
		Ek_temp = 0.5*atoms[n]->mass*dot(atoms[n]->v,atoms[n]->v);
		Ep_temp = atoms[n]->potential_energy;
		Ek += Ek_temp;
		Ep += Ep_temp;
		E += Ek_temp + Ep_temp;
	}

	fprintf(this->energyFile, "%f %f %f %f \n",t,Ek,Ep,E);
}

void StatisticsSampler::calculatePressure(double t) {
	if(!this->pressure) return;

	fprintf(this->pressureFile, "%f %f \n",t, this->system->P);
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

		fprintf(this->velocityFile, "%f %f %f %f \n",t, atom->v(0),atom->v(1),atom->v(2));
	}
}

void StatisticsSampler::calculateDiffusionConstant(double t) {
	
	if(!this->diffusionConstant || !t) return;

	int N = this->system->N;
	Atom **atoms = this->system->atoms;
	double rsquared = 0;
	for(int n=0;n<N;n++) {
		rsquared += atoms[n]->squaredDistanceFromInitialPosition();
	}
	rsquared /= N;
	double D = rsquared/(6*t);
	
	fprintf(this->diffusionFile, "%f %f\n",t, D);
	// fprintf(this->diffusionConstantFile, "%f %f\n",t, D);
	// fprintf(this->diffusionFile, "Nothing is happening");
}

