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

void Atom::update_velocity(const double &vx_, const double &vy_, const double &vz_) {
    // vx = vx_;
    // vy = vy_;
    // vz = vz_;
    v(0) = vx_;
    v(1) = vy_;
    v(2) = vz_;
}

void Atom::update_position(const double &rx_, const double &ry_, const double &rz_) {
    r(0) = rx_;
    r(1) = ry_;
    r(2) = rz_;
}

void Atom::update_initial_position(const double &rx_, const double &ry_, const double &rz_) {
    r_initial(0) = rx_;
    r_initial(1) = ry_;
    r_initial(2) = rz_;
}

void Atom::addR(const vec& dr) {
    r += dr;
    double L = system->L;

    for(int i=0;i<3;i++) {
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

const vec& Atom::distanceToAtom(Atom *atom, const vec &displacement) {
    dr = r-atom->r + displacement;
    return dr;
    /*
    double L = system->L;

	for(int i=0;i<3;i++) {
		if(dr(i) > L / 2.0) 
			dr(i) -= L;
		else if(dr(i)< -L / 2.0)
			dr(i) += L;
	}

	return dr;
    */
}

double Atom::squaredDistanceFromInitialPosition() {
    dr = r - r_initial;
    return dot(dr,dr);
}
