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
    double *r;
    double *v;
    double *r_initial;
    double *a;

    double vx, vy, vz, ax,ay,az, rx,ry,rz, r0x,r0y,r0z;

	System *system;

    Atom(System *system);
    ~Atom();
    void set_velocity(const double &vx_, const double &vy_, const double &vz_);
    void set_position(const double &rx_, const double &ry_, const double &rz_);
    void set_initial_position(const double &rx_, const double &ry_, const double &rz_);
    void step(const double &dt);
    void addR(const vec &dr);
    void distance_to_atom(Atom *atom, const vec &displacement, double &x, double &y, double &z);
	double squaredDistanceFromInitialPosition();
};
