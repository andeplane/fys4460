#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H

#include "Atom.h"
#include "System.h"
#include <fstream>

class StatisticsSampler {
private:
	System *system;
public:
	bool temperature;
	bool pressure;
	ofstream *temperatureFile;
	ofstream *pressureFile;
	ofstream *energyFile;

	StatisticsSampler(System *system);
	void sample(double t);
	void calculateTemperature(double t);
	void calculatePressure(double t);
	void calculateEnergy(double t);
};

#endif