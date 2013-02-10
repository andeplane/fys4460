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
    potential_energy = 0;
    dr = zeros<vec> (3,1);
}

Atom::~Atom() {

}

void Atom::addR(const vec& dr) {
    r += dr;
    double L = system->L;

    for(int i=0;i<3;i++) {
        int factor = r(i) / L;

        if(r(i) >= L) {
            int factor = r(i) / L;
            r(i) -= factor*L;
            r_initial(i) -= factor*L;
        }

        if(r(i) < 0) {
            int factor = r(i) / L;
            r(i) += (factor+1)*L;
            r_initial(i) = (factor+1)*L;
        }
    }
}

vec& Atom::distanceToAtom(Atom *atom) {
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
    dr = r - r_initial;
    return dot(dr,dr);
}
