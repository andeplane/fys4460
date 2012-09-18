#include "Atom.h"
#include <iostream>
#include <math.h>

using namespace std;

Atom::Atom(System *system) {
	this->r = zeros<vec> (3,1);
	this->v = zeros<vec> (3,1);
	this->F = zeros<vec> (3,1);
	this->mass = 39.948;         // MD units
	this->type = 0;
	this->system = system;
}

void Atom::addR(vec dr,double L) {
	this->r += dr;

	this->r(0) = fmod(this->r(0),L);
	this->r(1) = fmod(this->r(1),L);
	this->r(2) = fmod(this->r(2),L);
}

vec Atom::calculateForce(int startAt) {
	if(startAt==1) this->F.zeros(); // Reset force if this is a new calculation
	vec dr = zeros<vec>(3,1);
	vec F;
	Atom **atoms = this->system->getAtoms();
	Atom *atom;
	float drsquared;

	for(int n=startAt;n<this->system->getNumberOfAtoms();n++) {
		atom = atoms[n];
		if(startAt==1) atom->F.zeros(); // Reset force if this is a new calculation
		dr = this->r - atom->r;
		drsquared = norm(dr,1);

		F = 24*dr*(2/pow(drsquared,7) - 1/pow(drsquared,4));
		this->F += F;
		atom->F -= F;
	}

	return this->F;
}