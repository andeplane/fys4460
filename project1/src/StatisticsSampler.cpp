#include <statisticssampler.h>
#include <mdio.h>
#include <settings.h>
#include <unitconverter.h>
#include <system.h>
#include <mpi.h>
#include <mdtimer.h>
#include <iomanip.h>

StatisticsSampler::StatisticsSampler(System *system_) {
    system = system_;
    settings = system->settings;
    temperature_sampled_at = -1;
    kinetic_energy_sampled_at = -1;
    potential_energy_sampled_at = -1;
    pressure_sampled_at = -1;
}

void StatisticsSampler::sample_kinetic_energy() {
    if(system->steps == kinetic_energy_sampled_at) return;

    kinetic_energy = 0;
    double argon_mass = 39.948;
    double kinetic_energy_global = 0;

    for(unsigned int i=0;i<system->num_atoms_local;i++) {
        kinetic_energy += 0.5*argon_mass*(system->velocities[3*i+0]*system->velocities[3*i+0] + system->velocities[3*i+1]*system->velocities[3*i+1] + system->velocities[3*i+2]*system->velocities[3*i+2]);
    }
    MPI_Allreduce(&kinetic_energy, &kinetic_energy_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    kinetic_energy = kinetic_energy_global;
    kinetic_energy_sampled_at = system->steps;
}

void StatisticsSampler::sample_potential_energy() {
    if(system->steps == potential_energy_sampled_at) return;

    potential_energy = 0;
    MPI_Reduce(&system->potential_energy, &potential_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    potential_energy_sampled_at = system->steps;
}

void StatisticsSampler::sample_temperature() {
    if(system->steps == temperature_sampled_at) return;
    sample_kinetic_energy();
    double kinetic_energy_per_atom = kinetic_energy / system->num_atoms_global;
    temperature = 2.0/3*kinetic_energy_per_atom;

    temperature_sampled_at = system->steps;
}

void StatisticsSampler::sample_pressure() {
    if(system->steps == pressure_sampled_at) return;
    sample_temperature();

    pressure = 0;
    MPI_Reduce(&system->pressure_forces,&pressure,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    if(system->myid == 0) {
        pressure /= 3*system->volume;
        pressure += system->num_atoms_global/system->volume*temperature;
    }

    pressure_sampled_at = system->steps;
}

void StatisticsSampler::sample() {
    if(!settings->statistics_interval || system->steps % settings->statistics_interval != 0) return;

    system->mdtimer->start_sampling();
    double t_in_pico_seconds = system->unit_converter->time_to_SI(system->t)*1e12;

    sample_temperature();
    sample_potential_energy();
    sample_pressure();

    if(system->myid == 0) {
        double potential_energy_per_atom = potential_energy/system->num_atoms_global;
        double kinetic_energy_per_atom = kinetic_energy/system->num_atoms_global;

        fprintf(system->mdio->energy_file, "%.15f %.15f %.15f %.15f\n",t_in_pico_seconds,
                system->unit_converter->energy_to_ev(kinetic_energy_per_atom),
                system->unit_converter->energy_to_ev(potential_energy_per_atom),
                system->unit_converter->energy_to_ev(kinetic_energy_per_atom+potential_energy_per_atom),
                system->unit_converter->temperature_to_SI(temperature)
                );

        fprintf(system->mdio->pressure_file, "%.15f %.15f\n",t_in_pico_seconds,
                system->unit_converter->pressure_to_SI(pressure)
                );
        cout.setf(ios::fixed);
        cout.precision(5);
        cout << "Timestep " << setw(6) << system->steps << "   t=" << t_in_pico_seconds << " ps   T=" << system->unit_converter->temperature_to_SI(temperature) << " K" << endl;
    }

    system->mdtimer->end_sampling();
}
