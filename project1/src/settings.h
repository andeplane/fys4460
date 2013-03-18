#pragma once

#include <cinifile.h>

class Settings
{
public:
    Settings(string filename);
    CIniFile ini_file;
    double temperature;
    double dt;
    double FCC_b;
    double r_cut;
    int unit_cells_x;
    int unit_cells_y;
    int unit_cells_z;
    int nodes_x;
    int nodes_y;
    int nodes_z;
    int timesteps;
    int num_atoms_max;
    int movie_every_n_frame;
    bool create_movie;

};
