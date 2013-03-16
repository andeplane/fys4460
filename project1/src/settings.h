#pragma once

#include <CIniFile.h>

class Settings
{
public:
    Settings(string filename);
    CIniFile ini_file;
    double temperature;
    double dt;
    int unit_cells_x;
    int unit_cells_y;
    int unit_cells_z;
    int nodes_x;
    int nodes_y;
    int nodes_z;
    int timesteps;
    int max_particle_num;
};
