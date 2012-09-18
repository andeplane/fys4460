#ifndef ATOM_H
#define ATOM_H

#include <iostream>
#include <armadillo>
#include "System.h"

using namespace arma;
using namespace std;

class Atom {
public:
	char   type;
	double mass;
	vec r;
	vec v;
	vec a;
	vec F;
	System *system;
	

	Atom(System *system);
	void addR(vec dr,double L);
	vec calculateForce(int startAt);
};

#endif