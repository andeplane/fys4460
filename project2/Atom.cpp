#include "Atom.h"
#include <iostream>
#include <math.h>

using namespace std;

Atom::Atom(System *system) {
	this->r = zeros<vec> (3,1);
	this->v = zeros<vec> (3,1);
	this->F = zeros<vec> (3,1);
	this->a = zeros<vec> (3,1);
	this->mass = 1; // 39.948;         // MD units
	this->type = 0;
	this->system = system;

	this->interactingParticles = 0;
	this->interactingParticlesList = new int[system->N];
}

void Atom::addR(vec dr) {
	this->r += dr;
	double L = this->system->L;

	this->r(0) = fmod(this->r(0),L);
	this->r(1) = fmod(this->r(1),L);
	this->r(2) = fmod(this->r(2),L);
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