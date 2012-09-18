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
	System *system = new System();

	system->printPositionsToFile();

	return 0;
}