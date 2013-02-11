#include <System.h>
#ifdef MPI_ENABLED
#include <mpi.h>
#endif
// #define UNIFORMVELOCITY

double vmax = 8;


void System::initialize(double dt) {
    rnd = new Random(-rank - 1);
    N = 0;

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
        sort_cells();
    }

    if(rank==0) {
        for(int n=0;n<N;n++) {
            atoms[n]->v += 0.5*atoms[n]->a*dt;
            atoms[n]->addR(atoms[n]->v*dt);
        }
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
    r_and_v = new double[6*N];

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
