#pragma once

#include <iostream>
#include <armadillo>
#include "System.h"

using namespace arma;
using namespace std;

class Atom {
public:
	double mass;
    double *r;
    double *v;
    double *r_initial;
    double *a;
    int cell_index;
    Atom *prev, *next;
	System *system;

    Atom(System *system);
    ~Atom();
    void set_velocity(const double &vx, const double &vy, const double &vz);
    void set_position(const double &rx, const double &ry, const double &rz);
    void set_initial_position(const double &rx, const double &ry, const double &rz);
    void set_acceleration(const double &ax, const double &ay, const double &az);
    void step(const double &dt);
	double squaredDistanceFromInitialPosition();
};
