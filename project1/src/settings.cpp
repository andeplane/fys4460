#include "settings.h"

Settings::Settings(string filename) {
    ini_file.load(filename);

    temperature = ini_file.getdouble("temperature");
    dt = ini_file.getdouble("dt");
    unit_cells_x = ini_file.getint("unit_cells_x");
    unit_cells_y = ini_file.getint("unit_cells_y");
    unit_cells_z = ini_file.getint("unit_cells_z");
    nodes_x = ini_file.getint("nodes_x");
    nodes_y = ini_file.getint("nodes_y");
    nodes_z = ini_file.getint("nodes_z");
    timesteps = ini_file.getint("timesteps");
    max_particle_num = ini_file.getint("max_particle_num");
    print_positions = ini_file.getbool("print_positions");
}
