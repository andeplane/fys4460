#ifndef ATOM_H
#define ATOM_H

#include <iostream>
#include <armadillo>
#include "System.h"

using namespace arma;
using namespace std;

class Atom {
public:
	char   type;
	double mass;
	double potential_energy;
	vec r;
	vec v;
	vec a;
	vec F;
	System *system;
	int interactingParticles;
	int *interactingParticlesList;

	Atom(System *system);
	void addR(vec dr);
	vec distanceToAtom(Atom *atom);
};

#endif