#pragma once

class Atom;
class ThreadControl;
class Settings;
class MDIO;

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
    void calculate_accelerations();
    void half_kick();
    void kick();
    void move();
public:
    Settings *settings;
    MDIO *mdio;
    ThreadControl *thread_control;
    Random *rnd;
    int myid;
    double Lx, Ly, Lz;
    double dt;
    long steps;
    double Ek, Ep, E;

    System(int myid, Settings *settings);
    void step();
};
