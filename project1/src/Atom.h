#pragma once

#include <iostream>
#include <armadillo>
#include "System.h"

using namespace arma;
using namespace std;

class Atom {
public:
	double mass;
	double potential_energy;
    double *r;
    double *v;
    double *r_initial;
    double *a;
    int index_in_cell;
    int cell_index;

    double vx, vy, vz, ax,ay,az, rx,ry,rz, r0x,r0y,r0z;

	System *system;

    Atom(System *system);
    ~Atom();
    void set_velocity(const double &vx_, const double &vy_, const double &vz_);
    void set_position(const double &rx_, const double &ry_, const double &rz_);
    void set_initial_position(const double &rx_, const double &ry_, const double &rz_);
    void step(const double &dt);
	double squaredDistanceFromInitialPosition();
};
