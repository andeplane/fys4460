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
    void move();
public:
    vector<Atom*> atoms;
    vector<Atom*> all_atoms;
    Settings *settings;
    ThreadControl *thread_control;
    Random *rnd;
    int myid;
    double Lx, Ly, Lz;
    double dt;

    System(int myid, Settings *settings);
    void step();
};
