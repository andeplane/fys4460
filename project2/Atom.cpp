#include "Atom.h"
#include <iostream>
#include <math.h>

using namespace std;

Atom::Atom(System *system) {
	this->r = zeros<vec> (3,1);
	this->v = zeros<vec> (3,1);
	this->F = zeros<vec> (3,1);
	this->a = zeros<vec> (3,1);
	this->mass = 1; // 39.948;         // MD units
	this->type = 0;
	this->system = system;
}

void Atom::addR(vec dr,double L) {
	this->r += dr;

	this->r(0) = fmod(this->r(0),L);
	this->r(1) = fmod(this->r(1),L);
	this->r(2) = fmod(this->r(2),L);
}