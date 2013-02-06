#pragma once

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

	vec v;
	vec a;
	vec F;
	System *system;
	int interactingParticles;
	int *interactingParticlesList;
    Atom *next;
    Atom *prev;

	Atom(System *system);
    void addR(const vec dr);
	vec distanceToAtom(Atom *atom);
	double squaredDistanceFromInitialPosition();
};
