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
    vec temp_r;

    double vx, vy, vz, ax,ay,az, rx,ry,rz, r0x,r0y,r0z;

	vec v;
	vec a;
    vec dr;
    bool initialized;
    int index;
    int index2;

	System *system;
    Atom *next;
    Atom *prev;

    Atom(System *system);
    ~Atom();
    void update_velocity(const double &vx_, const double &vy_, const double &vz_);
    void update_position(const double &rx_, const double &ry_, const double &rz_);
    void update_initial_position(const double &rx_, const double &ry_, const double &rz_);
    void step(const double &dt);
    void addR(const vec &dr);
    void distance_to_atom(Atom *atom, const vec &displacement, double &x, double &y, double &z);
	double squaredDistanceFromInitialPosition();
};
