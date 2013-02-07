#include <System.h>
#include <mpi.h>
// #define UNIFORMVELOCITY

double vmax = 8;


void System::initialize() {
    rnd = new Random(-rank - 1);
    if(rank == 0) {
        init_atoms();
        printf("Initializing system with properties:\n");
        printf("T=%.2f\n",T);
        printf("L=%.2f\n",L);
        printf("Created %d atoms...",N);
    }

    // Send info about system size to all nodes
    MPI_Bcast(&L,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    V = L*L*L;

    cells_x = L/3;
    cells_y = L/3;
    cells_z = L/3;
    cell_width = L/cells_x;

    init_cells();

    if(rank==0) {
        send_particles_to_slaves();
    }
    else {
        receive_particles_from_master();
        //calculateAccelerations();
        cout << "I am " << rank << " and I calculated forces." << endl;
    }
}

void System::init_cells() {
    for(int k=0;k<cells_z;k++) {
        for(int j=0;j<cells_y;j++) {
            for(int i=0;i<cells_x;i++) {
                Cell *c = new Cell();
                c->i = i;
                c->j = j;
                c->k = k;
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
    vec zero_vec = zeros<vec>(3,1);
    double max_coord = 0;

	int n = 0;
    for(int x = 0; x < number_of_FCC_cells; x++) {
        for(int y = 0; y < number_of_FCC_cells; y++) {
            for(int z = 0; z < number_of_FCC_cells; z++) {
				for(int k = 0; k < 4; k++) {
                    Atom *atom = new Atom(this);

                    // Set positions and type
                    atom->r(0) = (x+xCell[k]) * b;
                    atom->r(1) = (y+yCell[k]) * b;
                    atom->r(2) = (z+zCell[k]) * b;
                    atom->r_initial = atom->r;

                    if(atom->r(0) > max_coord) max_coord = atom->r(0);
                    if(atom->r(1) > max_coord) max_coord = atom->r(1);
                    if(atom->r(2) > max_coord) max_coord = atom->r(2);

                    atom->type = k>0; // For visualization
                    atoms.push_back(atom);
                    atom->index = atoms.size()-1;
                }
            }
        }
    }

    L = b*number_of_FCC_cells;

    N = atoms.size();
    printf("done");
    initVelocities();
}

double System::gasdev() {
	static bool available = false;
	static double gset;
	double fac, rsq, v1, v2;
	if(!available) {
        do {
            v1 = 2.0*rnd->nextDouble() - 1.0;
            v2 = 2.0*rnd->nextDouble() - 1.0;
			rsq = v1*v1+v2*v2;
        } while(rsq>= 1.0 || rsq == 0.0);

		fac = sqrt(-2.0*log(rsq)/rsq);
		gset = v1*fac;
		available = true;
		return v2*fac;
	}
	else {
		available = false;
		return gset;
	}
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

    for(int n=0; n<N; n++ ) {
        atoms[n]->v(0) = rnd->nextGauss()*sqrt(k*T);
        atoms[n]->v(1) = rnd->nextGauss()*sqrt(k*T);
        atoms[n]->v(2) = rnd->nextGauss()*sqrt(k*T);
    }

    /*
    for(int n=0;n<N;n++)
		for(int i=0;i<3;i++) 
            atoms[n]->v(i) = gasdev();
    */
    double temp = 0;

    for(int n=0;n<N;n++) {
        temp += dot(atoms[n]->v,atoms[n]->v);
    }

    cout << "Temperature: " << temp/(3*N) << endl;

	vec vCM = zeros<vec>(3,1);
	
    for(int n=0;n<N;n++)
        vCM += atoms[n]->v;
		
    vCM /= N;
	
    for(int n=0;n<N;n++)
        atoms[n]->v -= vCM;
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
