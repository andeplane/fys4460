#ifndef ATOM_H
#define ATOM_H

#include <iostream>
using namespace std;

class Atom {
public:
	double *r;
	double *v;
	double *a;
	char   type;
	Atom();
	// string positionToString();
};

#endif