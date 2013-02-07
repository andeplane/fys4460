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
	double L;   // Length
    double V;
	double P;	// Pressure
    int cells_x, cells_y, cells_z;
    double cell_width;
    int rank;
    int nodes;
    ofstream *file;

    Random *rnd;
    vector<Cell*> cells;

	void printPositionsToFile(ofstream *file);
    void sort_cells();
    void send_particles_to_slaves();
    void receive_particles_from_master();
    void send_particles_back_to_master();
    void receive_particles_back_from_slaves();

	void step(double dt);
    System(int rank_, int nodes_, int number_of_FCC_cells_=4, double T_ = 1.0);

};