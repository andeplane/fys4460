#include "Atom.h"
#include <iostream>
#include <math.h>
#include <cstddef>

using namespace std;
int atom_count = 0;

Atom::Atom(System *system_) {
    initialized = true;
    next = NULL;
    prev = NULL;
    index2 = 0;
    r = zeros<vec> (3,1);
    v = zeros<vec> (3,1);
    a = zeros<vec> (3,1);

    mass = 1; // 39.948;         // MD units
    type = 0;
    system = system_;
}

Atom::~Atom() {

}

void Atom::addR(const vec& dr) {
    r += dr;
    double L = system->L;

    while(r(0) >= L) {
        r(0) -= L; r_initial(0) -= L;
    }
    while(r(0) < 0) {
        r(0) += L; r_initial(0) += L;
    }

    while(r(1) >= L)  {
        r(1) -= L; r_initial(1) -= L;
    }
    while(r(1) < 0)  {
        r(1) += L; r_initial(1) += L;
    }

    while(r(2) >= L)  {
        r(2) -= L; r_initial(2) -= L;
    }
    while(r(2) < 0)  {
        r(2) += L; r_initial(2) += L;
    }
}

vec Atom::distanceToAtom(Atom *atom) {
    dr = r-atom->r;

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
