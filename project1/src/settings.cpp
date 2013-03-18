#include <settings.h>
#include <mpi.h>
Settings::Settings(string filename) {
    try {
        ini_file.load(filename);

        FCC_b = ini_file.getdouble("FCC_b");
        temperature = ini_file.getdouble("temperature");
        dt = ini_file.getdouble("dt");
        r_cut = ini_file.getdouble("r_cut");
        unit_cells_x = ini_file.getint("unit_cells_x");
        unit_cells_y = ini_file.getint("unit_cells_y");
        unit_cells_z = ini_file.getint("unit_cells_z");
        nodes_x = ini_file.getint("nodes_x");
        nodes_y = ini_file.getint("nodes_y");
        nodes_z = ini_file.getint("nodes_z");
        timesteps = ini_file.getint("timesteps");
        num_atoms_max = ini_file.getint("num_atoms_max");
        movie_every_n_frame = ini_file.getint("movie_every_n_frame");
        create_movie = ini_file.getbool("create_movie");
    }
    catch (int e) {
        cout << "Could not load settings file." << endl;
        MPI_Finalize();
        exit(1);
    }
}
