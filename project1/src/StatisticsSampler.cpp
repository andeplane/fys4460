#include "StatisticsSampler.h"
#include <armadillo>

using namespace arma;

static int steps = 0;

StatisticsSampler::StatisticsSampler(System *system_) {
    system = system_;
    temperature = true;
    pressure = true;
    energy = true;
    printVelocities = false;
    diffusionConstant = true;

    temperatureFile = fopen("temperature.dat","w");
    pressureFile    = fopen("pressure.dat","w");
    energyFile      = fopen("energy.dat","w");
    velocityFile = fopen("velocity.dat","w");
    diffusionFile = fopen("diffusion.dat","w");
}

void StatisticsSampler::sample(double t) {
    calculateTemperature(t);
    calculateEnergy(t);
    calculatePressure(t);
    calculateVelocities(t);
    calculateDiffusionConstant(t);

	steps++;
}

void StatisticsSampler::calculateTemperature(double t) {
    if(!temperature) return;
    int N = system->N;
    double T = 0;

	for(int n=0;n<N;n++) {
        T += system->atoms[n]->v(0)*system->atoms[n]->v(0) + system->atoms[n]->v(1)*system->atoms[n]->v(1) + system->atoms[n]->v(2)*system->atoms[n]->v(2);
	}

    T/=3*N;

    fprintf(temperatureFile, "%f %f \n",t, T);
}

void StatisticsSampler::calculateEnergy(double t) {
    if(!energy) return;
    int N = system->N;
	double E = 0, Ek=0,Ep=0, Ek_temp, Ep_temp;

    for(int n=0;n<N;n++) {
        Ek_temp = 0.5*system->atoms[n]->mass*dot(system->atoms[n]->v,system->atoms[n]->v);
        Ep_temp = system->atoms[n]->potential_energy;
		Ek += Ek_temp;
		Ep += Ep_temp;
		E += Ek_temp + Ep_temp;
	}

    // cout << "energy = " << E << endl;

    fprintf(energyFile, "%f %f %f %f \n",t,Ek,Ep,E);
}

void StatisticsSampler::calculatePressure(double t) {
    if(!pressure) return;

    fprintf(pressureFile, "%f %f \n",t, system->P);
}

void StatisticsSampler::calculateVelocities(double t) {
    if(!printVelocities || steps % 50) return;

    int N = system->N;

    if(!steps) fprintf(velocityFile, "%d %d %d\n",N,N,N);
	Atom *atom;
	double vsq = 0;

	for(int n=0;n<N;n++) {
        atom = system->atoms[n];
		// vsq = norm(atom->v,2);

        fprintf(velocityFile, "%f %f %f %f \n",t, atom->v(0),atom->v(1),atom->v(2));
	}
}

void StatisticsSampler::calculateDiffusionConstant(double t) {
	
    if(!diffusionConstant || !t) return;

    int N = system->N;

	double rsquared = 0;
	for(int n=0;n<N;n++) {
        rsquared += system->atoms[n]->squaredDistanceFromInitialPosition();
	}
	rsquared /= N;
	double D = rsquared/(6*t);
	
    fprintf(diffusionFile, "%f %f\n",t, D);
    // fprintf(diffusionConstantFile, "%f %f\n",t, D);
    // fprintf(diffusionFile, "Nothing is happening");
}

