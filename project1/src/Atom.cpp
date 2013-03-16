#include "Atom.h"
#include <iostream>
#include <math.h>
#include <cstddef>

using namespace std;
int atom_count = 0;

Atom::Atom(System *system_) {
    r_initial = NULL;
    r = NULL;
    v = NULL;
    a = NULL;

    mass = 1; // 39.948;         // MD units
    system = system_;
    potential_energy = 0;
}

Atom::~Atom() {

}

void Atom::set_velocity(const double &vx, const double &vy, const double &vz) {
    v[0] = vx;
    v[1] = vy;
    v[2] = vz;
}

void Atom::set_position(const double &rx, const double &ry, const double &rz) {
    r[0] = rx;
    r[1] = ry;
    r[2] = rz;
}

void Atom::set_initial_position(const double &rx, const double &ry, const double &rz) {
    r_initial[0] = rx;
    r_initial[1] = ry;
    r_initial[2] = rz;
}

void Atom::step(const double &dt) {
    r[0] += v[0]*dt;
    r[1] += v[1]*dt;
    r[2] += v[2]*dt;

    double Lx = system->Lx;
    double Ly = system->Ly;
    double Lz = system->Lz;

    if(r[0]>=Lx) {
        int factor = r[0] / Lx;
        r[0] -= factor*Lx;
        r_initial[0] -= factor*Lx;
    } else if(r[0] < 0) {
        int factor = r[0] / Lx;
        r[0] += (factor+1)*Lx;
        r_initial[0] = (factor+1)*Lx;
    }

    if(r[1]>=Ly) {
        int factor = r[1] / Ly;
        r[1] -= factor*Ly;
        r_initial[1] -= factor*Ly;
    } else if(r[1] < 0) {
        int factor = r[1] / Ly;
        r[1] += (factor+1)*Ly;
        r_initial[1] = (factor+1)*Ly;
    }

    if(r[2]>=Lz) {
        int factor = r[2] / Lz;
        r[2] -= factor*Lz;
        r_initial[2] -= factor*Lz;
    } else if(r[2] < 0) {
        int factor = r[2] / Lz;
        r[2] += (factor+1)*Lz;
        r_initial[2] = (factor+1)*Lz;
    }
}

double Atom::squaredDistanceFromInitialPosition() {
    return pow(r[0] - r_initial[0],2) + pow(r[1] - r_initial[1],2) + pow(r[2] - r_initial[2],2);
}
