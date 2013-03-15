#pragma once

class Atom;
class ThreadControl;
class Settings;

#include <fstream>
#include <Atom.h>
#include <random.h>
#include <threadcontrol.h>
#include <Cell.h>
#include <vector>

using namespace std;

class System {
private:
    void initialize();
	void calculateAccelerations();
    void update_velocity_and_move(const double &dt);
public:
    vector<Atom*> atoms;
    vector<Atom*> all_atoms;
    Settings *settings;
    ThreadControl *thread_control;
    Random *rnd;
    int myid;
    double Lx, Ly, Lz;

    void step();
    System(int myid, Settings *settings);
};
