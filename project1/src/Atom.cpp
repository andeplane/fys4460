#include "Atom.h"
#include <iostream>
#include <math.h>
#include <cstddef>

using namespace std;

Atom::Atom(System *system_) {
    next = NULL;
    prev = NULL;
    r = zeros<vec> (3,1);
    v = zeros<vec> (3,1);
    F = zeros<vec> (3,1);
    a = zeros<vec> (3,1);

    mass = 1; // 39.948;         // MD units
    type = 0;
    system = system_;

    interactingParticles = 0;
    interactingParticlesList = new int[system->N];
}

void Atom::addR(const vec dr) {
    r += dr;
    double L = system->L;

    if(r(0) > L)  { r(0) -= L; r_initial(0) -= L; }
    else if(r(0) < 0)  { r(0) += L; r_initial(0) += L; }

    if(r(1) > L)  { r(1) -= L; r_initial(1) -= L; }
    else if(r(1) < 0)  { r(1) += L; r_initial(1) += L; }

    if(r(2) > L)  { r(2) -= L; r_initial(2) -= L; }
    else if(r(2) < 0)  { r(2) += L; r_initial(2) += L; }
}

vec Atom::distanceToAtom(Atom *atom) {
    vec dr = r-atom->r;
    double L = system->L;

	for(int i=0;i<3;i++) {
		if(dr(i) > L / 2.0) 
			dr(i) -= L;
		else if(dr(i)< -L / 2.0)
			dr(i) += L;
	}

	return dr;
}

double Atom::squaredDistanceFromInitialPosition() {
    vec dr = r - r_initial;
	return dot(dr,dr);
}
