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
    vec dr;
    bool initialized;
    int index;
    int index2;

	System *system;
    Atom *next;
    Atom *prev;

	Atom(System *system);
    void addR(const vec &dr);
	vec distanceToAtom(Atom *atom);
	double squaredDistanceFromInitialPosition();
};
