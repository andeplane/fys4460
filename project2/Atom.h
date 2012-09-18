#ifndef ATOM_H
#define ATOM_H

#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;

class Atom {
public:
	double mass;
	vec r;
	vec v;
	vec a;
	
	char   type;
	Atom();
};

#endif