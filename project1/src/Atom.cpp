#include "Atom.h"
#include <iostream>
#include <math.h>

using namespace std;

Atom::Atom(System *system_) {
    r = zeros<vec> (3,1);
    v = zeros<vec> (3,1);
    F = zeros<vec> (3,1);
    a = zeros<vec> (3,1);
    temp_vector = zeros<vec> (3,1);
    mass = 1; // 39.948;         // MD units
    type = 0;
    system = system_;

    interactingParticles = 0;
    interactingParticlesList = new int[system->N];
}

void Atom::addR(vec dr) {
    r += dr;
    double L = system->L;

    temp_vector(0) = fmod(r(0)+10*L,L);
    temp_vector(1) = fmod(r(1)+10*L,L);
    temp_vector(2) = fmod(r(2)+10*L,L);
	
    r_initial(0) -= (r(0) - temp_vector(0));
    r_initial(1) -= (r(1) - temp_vector(1));
    r_initial(2) -= (r(2) - temp_vector(2));

    r = temp_vector;
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
