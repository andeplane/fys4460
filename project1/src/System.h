#pragma once

class Atom;
class Cell;

#include <fstream>
#include "Atom.h"
#include <random.h>
#include <Cell.h>

using namespace std;

class System {
private:
	void initialize();
    void init_cells();
	void initPositions();
	void initVelocities();
	void updateVerletList();
	void rescaleVelocities();
	void calculateAccelerations();
	
	double gasdev();

public:
	Atom **atoms;
	int N; 		// Number of atoms
	double T; 	// Temperature
	double rho; // Density
	double L;   // Length
	double P;	// Pressure
    int cells_x, cells_y, cells_z;

    Random *rnd;
    vector<Cell> cells;

	void printPositionsToFile(ofstream *file);
    void sort_cells();

	void step(double dt);
	System(int N=108, double T=1.0, double rho=0.8);

};
