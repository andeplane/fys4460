#pragma once

#include <stdio.h>
#include <system.h>
#include <fstream>

class System;
class Settings;

class StatisticsSampler {
private:
    System *system;
    Settings *settings;
public:
    StatisticsSampler(System *system);
    void sample();
    double kinetic_energy, potential_energy, temperature, mean_r_squared;
};
