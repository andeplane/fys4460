/*
#include "thermostat.h"
#include <vector>
using namespace std;

Thermostat::Thermostat(double T_, double relaxation_time_, double dt_)
{
    T = T_;
    relaxation_time = relaxation_time_;
    dt = dt_;
}

void Thermostat::apply(const vector<Atom*>& atoms) {
    double system_temperature = 0;
    for(int n=0;n<atoms.size();n++) {
        system_temperature += dot(atoms[n]->v,atoms[n]->v);
    }
    system_temperature /= 3*atoms.size();

    double berendsen_factor = sqrt(1 + dt/relaxation_time*(T/system_temperature - 1));

    for(int n=0;n<atoms.size();n++) {
        atoms[n]->v *= berendsen_factor;
    }
}
*/
