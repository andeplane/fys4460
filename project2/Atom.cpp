#include "Atom.h"
#include <iostream>
#include <math.h>

using namespace std;

Atom::Atom() {
	this->r = zeros<vec> (3,1);
	this->v = zeros<vec> (3,1);
	this->a = zeros<vec> (3,1);
	this->mass = 39.948;         // MD units
	this->type = 0;
}