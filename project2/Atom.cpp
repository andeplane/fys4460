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

void Atom::addR(vec dr,double L) {
	this->r += dr;

	this->r(0) = fmod(this->r(0),L);
	this->r(1) = fmod(this->r(1),L);
	this->r(2) = fmod(this->r(2),L);
}