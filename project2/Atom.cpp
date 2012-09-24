#include "Atom.h"
#include <iostream>
#include <math.h>

using namespace std;

Atom::Atom(System *system) {
	this->r = zeros<vec> (3,1);
	this->v = zeros<vec> (3,1);
	this->F = zeros<vec> (3,1);
	this->a = zeros<vec> (3,1);
	this->temp_vector = zeros<vec> (3,1);
	this->mass = 1; // 39.948;         // MD units
	this->type = 0;
	this->system = system;

	this->interactingParticles = 0;
	this->interactingParticlesList = new int[system->N];
}

void Atom::addR(vec dr) {
	this->r += dr;
	double L = this->system->L;

	this->temp_vector(0) = fmod(this->r(0)+10*L,L);
	this->temp_vector(1) = fmod(this->r(1)+10*L,L);
	this->temp_vector(2) = fmod(this->r(2)+10*L,L);
	
	this->r_initial(0) -= (this->r(0) - this->temp_vector(0));
	this->r_initial(1) -= (this->r(1) - this->temp_vector(1));
	this->r_initial(2) -= (this->r(2) - this->temp_vector(2));

	this->r = this->temp_vector;
}

vec Atom::distanceToAtom(Atom *atom) {
	vec dr = this->r-atom->r;
	double L = this->system->L;

	for(int i=0;i<3;i++) {
		if(dr(i) > L / 2.0) 
			dr(i) -= L;
		else if(dr(i)< -L / 2.0)
			dr(i) += L;
	}

	return dr;
}

double Atom::squaredDistanceFromInitialPosition() {
	vec dr = this->r - this->r_initial;
	return dot(dr,dr);
}