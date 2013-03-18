#include <statisticssampler.h>
#include <mdio.h>
#include <settings.h>
#include <unitconverter.h>
#include <system.h>
#include <mpi.h>
#include <mdtimer.h>

StatisticsSampler::StatisticsSampler(System *system_) {
    system = system_;
    settings = system->settings;

}

void StatisticsSampler::sample() {
    system->mdtimer->start_sampling();

    kinetic_energy = 0;
    double argon_mass = 39.948;

    if(settings->statistics_interval && system->steps % settings->statistics_interval != 0) return;
    for(unsigned int i=0;i<system->num_atoms_local;i++) {
        kinetic_energy += 0.5*argon_mass*system->velocities[3*i+0]*system->velocities[3*i+0] + system->velocities[3*i+1]*system->velocities[3*i+1] + system->velocities[3*i+2]*system->velocities[3*i+2];
    }

    double t_in_pico_seconds = system->unit_converter->time_to_SI(system->t)*1e12;
    double kinetic_energy_global = 0;
    potential_energy = 0;

    MPI_Reduce(&kinetic_energy, &kinetic_energy_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&system->potential_energy, &potential_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    kinetic_energy = kinetic_energy_global; // Update local variable

    if(system->myid == 0) {
        temperature = 2.0/3*kinetic_energy/system->num_atoms_global;

        fprintf(system->mdio->energy_file, "%f %f %f %f\n",t_in_pico_seconds,
                system->unit_converter->energy_to_ev(kinetic_energy),
                system->unit_converter->energy_to_ev(potential_energy),
                system->unit_converter->energy_to_ev(kinetic_energy+potential_energy),
                system->unit_converter->temperature_to_SI(temperature)
                );
        cout << system->steps << "   t=" << t_in_pico_seconds << "   T=" << system->unit_converter->temperature_to_SI(temperature) << endl;
    }

    system->mdtimer->end_sampling();
}
