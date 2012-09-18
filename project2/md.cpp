#include <iostream>
#include "math.h"
#include "time.h"
#include <armadillo>
#include <fstream>
#include "Atom.h"
#include "System.h"

using namespace arma;
using namespace std;

int main(int argc, char *argv[]) {
	System *system = new System(256);

	ofstream *file = new ofstream;
	file->open("pos.xyz");

	system->printPositionsToFile(file);
	for(int i=0;i<500;i++) {
		system->step(0.005);
		system->printPositionsToFile(file);
	}
	
	file->close();

	return 0;
}