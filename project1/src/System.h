#pragma once

class Atom;
class Cell;

#include <fstream>
#include <Atom.h>
#include <random.h>
#include <Cell.h>
#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

class System {
private:
	void initialize();
    void init_cells();
    void init_atoms();
	void initVelocities();
	void rescaleVelocities();
	void calculateAccelerations();
	
	double gasdev();

public:
    vector<Atom*> atoms;
	int N; 		// Number of atoms
    int number_of_FCC_cells;
	double T; 	// Temperature
	double rho; // Density
	double L;   // Length
	double P;	// Pressure
    int cells_x, cells_y, cells_z;
    double cell_width;

    Random *rnd;
    vector<Cell*> cells;

	void printPositionsToFile(ofstream *file);
    void sort_cells();

	void step(double dt);
    System(int number_of_FCC_cells_=4, double T=1.0, double rho=0.8);

};
