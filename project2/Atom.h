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
	vec r_initial; // Initial position of the atom. This will be changed while going through the boundary
	vec temp_vector; // Used to temporarily save vectors
	vec v;
	vec a;
	vec F;
	System *system;
	int interactingParticles;
	int *interactingParticlesList;

	Atom(System *system);
	void addR(vec dr);
	vec distanceToAtom(Atom *atom);
	double squaredDistanceFromInitialPosition();
};

#endif