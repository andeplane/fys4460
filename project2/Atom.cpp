#include "Atom.h"
#include <iostream>
using namespace std;

Atom::Atom() {
	this->r = new double[3];
	this->v = new double[3];
	this->a = new double[3];
	this->type = 0;
}

// string Atom::positionToString() {
	// file << (type[n] ? "Ar " : "H ") << r[n][0] << " " << r[n][1] << " " << r[n][2] << endl;
	// return sprintf("%s %f %f %f",(this->type ? "H" : "Ar"),this->r[0],this->r[1],this->r[2]);
// }