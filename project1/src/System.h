#pragma once

class Atom;

#include <fstream>
#include "Atom.h"

using namespace std;

class System {
private:
	void initialize();
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
	
	// int getNumberOfAtoms();
	// double getTemperature();
	// double getDensity();
	// double getLength();
	void printPositionsToFile(ofstream *file);

	void step(double dt);
	System(int N=108, double T=1.0, double rho=0.8);

};
