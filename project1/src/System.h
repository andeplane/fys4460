#pragma once

class Atom;
class ThreadControl;
class Settings;
class MDIO;

#include <fstream>
#include <random.h>
#include <vector>

using namespace std;

class System {
private:
    void initialize();
    void calculate_accelerations();
    void half_kick();
    void kick();
    void atom_move();
    void atom_copy();
    void create_FCC();
public:
    Settings *settings;
    MDIO *mdio;
    Random *rnd;

    long steps;
    int myid;
    int node_idx, node_idy, node_idz;

    int num_nodes;
    unsigned long num_particles_local;
    unsigned long num_particles_global;
    unsigned long num_particles_ghost;
    double r_cut, Lx, Ly, Lz, dt;
    double Lx_full, Ly_full, Lz_full;
    double *positions;
    double *accelerations;
    double *velocities;
    double origo[3];

    System();
    void setup(int myid_, Settings *settings_);
    void step();
};
