#include <System.h>
#include <unitconverter.h>
#include <settings.h>

#include <mpi.h>
// #define UNIFORMVELOCITY

double vmax = 8;

void calc_stuff(System *system) {
    double dr_2, dr_6, dr_12, f, potential_energy, dr_12_inv, dr_6_inv;
    double x,y,z;
    double L = system->L;
    int N = system->N;
    Atom *atom0, *atom1;
    for(int i=0;i<N;i++) {
        atom0 = system->atoms[i];
        for(int j=i+1;j<N;j++) {
            atom1 = system->atoms[j];
            // atom0->distance_to_atom(atom1,displacement_vector,x,y,z);
            x = atom0->r(0) - atom1->r(0);
            y = atom0->r(1) - atom1->r(1);
            z = atom0->r(2) - atom1->r(2);

            if(x>L/2) {
                x -= L;
            }

            if(y>L/2) {
                y -= L;
            }

            if(z>L/2) {
                z -= L;
            }

            dr_2 = x*x + y*y + z*z;
            if(dr_2>9) continue;

            dr_6 = pow(dr_2,3);
            dr_12 = pow(dr_6,2);
            dr_12_inv = 1.0/dr_12;
            dr_6_inv = 1.0/dr_6;


            f = 24*(2.0*dr_12_inv-dr_6_inv)/dr_2;

            potential_energy = 4*(dr_12_inv - dr_6_inv);

            atom0->potential_energy += potential_energy;
        }
    }
}

void System::initialize() {
    rnd = new Random(-(myid+1));
    double L = 1.54478707783;

    Lx = settings->unit_cells_x*L;
    Ly = settings->unit_cells_y*L;
    Lz = settings->unit_cells_z*L;

    if(rank == 0) {
        init_atoms();
    }

#ifdef MPI_ENABLED
    // Send info about system size to all nodes
    MPI_Bcast(&L,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
    V = L*L*L;

    cells_x = L/3;
    cells_y = L/3;
    cells_z = L/3;
    cell_width = L/cells_x;

    init_cells();
    thread_control = new ThreadControl();

    thread_control->setup(nodes,cells_x,cells_y,cells_z,cells);

    if(rank==0) {
        printf("Initializing system...\n");
        printf("T=%.2f\n",T);
        printf("L=%.2f\n",L);
        printf("%d atoms",N);
        printf("%d cells",(int)cells.size());
        printf("%d cells per dim",cells_z);
        // atoms.clear();
    }

    UnitConverter uc;
    // calculateAccelerations();
    calc_stuff(this);

    double E = 0;
    for(int n=0;n<N;n++) {
        E += atoms[n]->potential_energy;
    }
    cout << "Potential energy: " << E << endl;
    E /= N;
    double E_SI = uc.energy_to_SI(E);

    cout << "energy per particle: " << E << endl;
    cout << "energy per particle (SI): " << E_SI << endl;

    exit(0);

    if(rank==0) {
        for(int n=0;n<N;n++) {
            atoms[n]->v += 0.5*atoms[n]->a*dt;
            // atoms[n]->addR(atoms[n]->v*dt);
        }
    }

    calculateAccelerations();
}

void System::init_cells() {
    for(int k=0;k<cells_z;k++) {
        for(int j=0;j<cells_y;j++) {
            for(int i=0;i<cells_x;i++) {
                Cell *c = new Cell();
                c->i = i;
                c->j = j;
                c->k = k;
                c->index = cells.size();
                cells.push_back(c);
            }
        }
    }

    for(int i=0;i<cells.size();i++) {
        cells[i]->find_neighbours(cells_x,cells_y,cells_z, this);
    }
}

void System::init_atoms() {
	printf("Initializing FCC lattice...");

    double b = 1.545;

    double xCell[4] = {0, 0.5, 0.5, 0};
    double yCell[4] = {0, 0.5, 0, 0.5};
    double zCell[4] = {0, 0, 0.5, 0.5};
    double rx,ry,rz;

    for(int x = 0; x < number_of_FCC_cells; x++) {
        for(int y = 0; y < number_of_FCC_cells; y++) {
            for(int z = 0; z < number_of_FCC_cells; z++) {
				for(int k = 0; k < 4; k++) {
                    Atom *atom = new Atom(this);

                    rx = (x+xCell[k]) * b;
                    ry = (y+yCell[k]) * b;
                    rz = (z+zCell[k]) * b;

                    atom->update_position(rx,ry,rz);
                    atom->update_initial_position(rx,ry,rz);
                    // Set positions and type

                    atom->type = k>0; // For visualization
                    atom->index = atoms.size();

                    atoms.push_back(atom);
                    all_atoms.push_back(atom);
                }
            }
        }
    }

    L = b*number_of_FCC_cells;

    N = atoms.size();
    printf("done");
    initVelocities();
}

void System::initVelocities() {
#ifdef UNIFORMVELOCITY
	printf("Creating uniform velocities...");
	double realVmax = vmax/sqrt(3);

    for(int n=0;n<N;n++) {
		for(int k=0;k<3;k++) 
            atoms[n]->v(k) = 2*ran0(&idum)-1;
        atoms[n]->v *= realVmax;
	}
	printf("done\n");
#else

	printf("Creating maxwellian velocities...");

    double k = 1.0;
    double vx,vy,vz;

    for(int n=0; n<N; n++ ) {
        vx = rnd->nextGauss()*sqrt(k*T);
        vy = rnd->nextGauss()*sqrt(k*T);
        vz = rnd->nextGauss()*sqrt(k*T);

        atoms[n]->update_velocity( vx , vy, vz );
    }

    double temp = 0;

    for(int n=0;n<N;n++) {
        temp += dot(atoms[n]->v,atoms[n]->v);
    }

    cout << "Temperature: " << temp/(3*N) << endl;

	vec vCM = zeros<vec>(3,1);
	
    for(int n=0;n<N;n++)
        vCM += atoms[n]->v;
		
    vCM /= N;
	
    for(int n=0;n<N;n++) {
        /*
        vx = r_and_v[6*n+3] - vCM(0);
        vy = r_and_v[6*n+4] - vCM(1);
        vz = r_and_v[6*n+5] - vCM(2);
        */
        // atoms[n]->update_velocity( vx , vy, vz );

        atoms[n]->v -= vCM;
    }

	printf("done\n");
    rescaleVelocities();
#endif
}

void System::rescaleVelocities() {
	printf("Rescaling velocities...");
	double vSqdSum = 0;
    for(int n=0;n<N;n++)
        vSqdSum += norm(atoms[n]->v,1);
    double lambda = sqrt(3*(N-1)*T/vSqdSum);
    for(int n=0;n<N;n++)
        atoms[n]->v *= lambda;
	printf("done\n");
}
