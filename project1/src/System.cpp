#include <iostream>
#include "math.h"
#include "time.h"
#include <fstream>
#include "System.h"
#include <inlines.h>
#include <settings.h>
#include <mpi.h>
#include <mdio.h>
#include <unitconverter.h>

using namespace std;
#define MOVED_OUT -1.0e10

System::System() {

}

void System::setup(int myid_, Settings *settings_) {
    myid = myid_;
    settings = settings_;
    dt = settings->dt;
    dt_half = dt/2;
    r_cut = settings->r_cut;
    num_atoms_local = 0;
    num_atoms_global = 0;
    num_atoms_ghost = 0;
    num_nodes = settings->nodes_x*settings->nodes_y*settings->nodes_z;

    mpi_send_buffer = new double[settings->num_atoms_max];
    mpi_receive_buffer = new double[settings->num_atoms_max];

    positions = new double[3*settings->num_atoms_max];
    accelerations = new double[3*settings->num_atoms_max];
    velocities = new double[3*settings->num_atoms_max];
    atom_moved = new bool[settings->num_atoms_max];
    for(i=0;i<6;i++) move_queue[i] = new int[settings->num_atoms_max];

    steps = 0;
    rnd = new Random(-(myid+1));

    num_processors[0] = settings->nodes_x;
    num_processors[1] = settings->nodes_y;
    num_processors[2] = settings->nodes_z;

    node_index[0] = myid/(settings->nodes_y*settings->nodes_z);
    node_index[1] = (myid/settings->nodes_z) % settings->nodes_y;
    node_index[2] = myid%settings->nodes_z;

    // Size of this node
    box_length[0] = settings->unit_cells_x*settings->FCC_b;
    box_length[1] = settings->unit_cells_y*settings->FCC_b;
    box_length[2] = settings->unit_cells_z*settings->FCC_b;

    for(a=0;a<3;a++) {
        box_length_full[a] = box_length[a]*num_processors[a];
        origo[a] = (float)node_index[a] * box_length[a];
    }

    mdio = new MDIO();
    mdio->setup(this);
    set_topology();
    create_FCC();
    mpi_copy();
    calculate_accelerations();
    half_kick();

    if(myid==0) cout << "System size: " << box_length_full[0] << " " << box_length_full[1] << " " << box_length_full[2] << endl;
    if(myid==0) cout << "Atoms: " << num_atoms_global << endl;
}

void System::create_FCC() {
    UnitConverter *unit_converter = new UnitConverter();
    double xCell[4] = {0, 0.5, 0.5, 0};
    double yCell[4] = {0, 0.5, 0, 0.5};
    double zCell[4] = {0, 0, 0.5, 0.5};

    double r[3];
    double T = unit_converter->temperature_from_SI(settings->temperature);

    int n = 0;
    for(int x = 0; x < settings->nodes_x*settings->unit_cells_x; x++) {
        for(int y = 0; y < settings->nodes_y*settings->unit_cells_y; y++) {
            for(int z = 0; z < settings->nodes_z*settings->unit_cells_z; z++) {
                for(int k = 0; k < 4; k++) {

                    // Set positions and type
                    r[0] = (x+xCell[k]) * settings->FCC_b - origo[0];
                    r[1] = (y+yCell[k]) * settings->FCC_b - origo[1];
                    r[2] = (z+zCell[k]) * settings->FCC_b - origo[2];
                    bool is_mine = true;
                    for(i=0;i<3;i++) {
                        if(!(r[i] >= 0 && r[i] < box_length[i])) is_mine = false;
                    }

                    if(is_mine) {
                        for(i=0;i<3;i++) positions[3*num_atoms_local+i] = r[i];

                        velocities[3*num_atoms_local+0] = rnd->nextGauss()*sqrt(T);
                        velocities[3*num_atoms_local+1] = rnd->nextGauss()*sqrt(T);
                        velocities[3*num_atoms_local+2] = rnd->nextGauss()*sqrt(T);

                        num_atoms_local++;
                    }
                }
            }
        }
    }

    MPI_Allreduce(&num_atoms_local,&num_atoms_global,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
}

void System::mpi_copy() {
    MPI_Status status;

    int node_id, num_send, num_receive;
    int new_ghost_atoms = 0;
    short higher, local_node_id;
    for(short dimension=0;dimension<3;dimension++) {
        for (higher=0; higher<2; higher++) move_queue[2*dimension+higher][0] = 0;
        for(i=0;i<num_atoms_local+new_ghost_atoms;i++) {
            for(higher=0;higher<2;higher++) {
                local_node_id = 2*dimension + higher;
                if (atom_should_be_copied(&positions[3*i],local_node_id)) move_queue[local_node_id][++(move_queue[local_node_id][0])] = i;
            }
        }

        /* Loop through higher and lower node in this dimension */
        for(higher=0;higher<2;higher++) {
            local_node_id= 2*dimension+higher;
            node_id = neighbor_nodes[local_node_id];
            num_send = move_queue[local_node_id][0];

            if (my_parity[dimension] == 0) {
                MPI_Send(&num_send,1,MPI_INT,node_id,10,MPI_COMM_WORLD);
                MPI_Recv(&num_receive,1,MPI_INT,MPI_ANY_SOURCE,10,MPI_COMM_WORLD,&status);
            }
            else if (my_parity[dimension] == 1) {
                MPI_Recv(&num_receive,1,MPI_INT,MPI_ANY_SOURCE,10,MPI_COMM_WORLD,&status);
                MPI_Send(&num_send,1,MPI_INT,node_id,10,MPI_COMM_WORLD);
            }
            else num_receive = num_send;

            for (i=1; i<=num_send; i++) {
                for (a=0; a<3; a++) { /* Shift the coordinate origin */
                  mpi_send_buffer[3*(i-1)+a] = positions[ 3*move_queue[local_node_id][i] + a]-shift_vector[local_node_id][a];
                }
            }

            if (my_parity[dimension] == 0) {
                MPI_Send(mpi_send_buffer,3*num_send,MPI_DOUBLE,node_id,20,MPI_COMM_WORLD);
                MPI_Recv(mpi_receive_buffer,3*num_receive,MPI_DOUBLE,MPI_ANY_SOURCE,20,MPI_COMM_WORLD,&status);
            }
            else if (my_parity[dimension] == 1) {
                MPI_Recv(mpi_receive_buffer,3*num_receive,MPI_DOUBLE,MPI_ANY_SOURCE,20,MPI_COMM_WORLD,&status);
                MPI_Send(mpi_send_buffer,3*num_send,MPI_DOUBLE,node_id,20,MPI_COMM_WORLD);
            }
            else for (i=0; i<3*num_receive; i++) mpi_receive_buffer[i] = mpi_send_buffer[i];

            for (i=0; i<num_receive; i++) {
                for (a=0; a<3; a++) positions[ 3*(num_atoms_local+new_ghost_atoms+i) + a] = mpi_receive_buffer[3*i+a];
            }

            new_ghost_atoms += num_receive;
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    num_atoms_ghost = new_ghost_atoms;
}

void System::set_topology() {
    /*----------------------------------------------------------------------
    Defines a logical network topology.  Prepares a neighbor-node ID table,
    nn, & a shift-vector table, sv, for internode message passing.  Also
    prepares the node parity table, myparity.
    ----------------------------------------------------------------------*/

    /* Integer vectors to specify the six neighbor nodes */
    int iv[6][3] = {
        {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}
    };

    int k1[3];

    /* Set up neighbor tables, nn & sv */
    for (int ku=0; ku<6; ku++) {
        /* Vector index of neighbor ku */
        for (a=0; a<3; a++) {
            k1[a] = (node_index[a]+iv[ku][a]+num_processors[a])%num_processors[a];
        }

        /* Scalar neighbor ID, nn */
        neighbor_nodes[ku] = k1[0]*num_processors[1]*num_processors[2]+k1[1]*num_processors[2]+k1[2];
        /* Shift vector, sv */
        for (a=0; a<3; a++) shift_vector[ku][a] = box_length[a]*iv[ku][a];
    }



    /* Set up the node parity table, myparity */
    for (a=0; a<3; a++) {
        if (num_processors[a] == 1) {
            my_parity[a] = 2;
        } else if (node_index[a]%2 == 0) {
            my_parity[a] = 0;
        } else {
            my_parity[a] = 1;
        }
    }
}

inline bool System::atom_should_be_copied(double* ri, int ku) {
  int dimension,higher;
  dimension = ku/2; /* x(0)|y(1)|z(2) direction */
  higher = ku%2; /* Lower(0)|higher(1) direction */
  if (higher == 0) return ri[dimension] < r_cut;
  else return ri[dimension] > box_length[dimension]-r_cut;
}


inline bool System::atom_did_change_node(double* ri, int ku) {
    int dimension,higher;
    dimension = ku/2;    /* x(0)|y(1)|z(2) direction */
    higher = ku%2; /* Lower(0)|higher(1) direction */
    if (higher == 0) return ri[dimension] < 0.0;
    else return ri[dimension] > box_length[dimension];
}

void System::mpi_move() {
    MPI_Status status;

    int new_atoms = 0;
    int node_id,num_send,num_receive;
    short node_higher, node_lower, local_node_id;

    double com1;

    /* Reset the # of to-be-moved atoms, move_queue[][0] */
    for (short ku=0; ku<6; ku++) move_queue[ku][0] = 0;

    /* Make a moved-atom list, move_queue----------------------------------*/

    for(short dimension=0;dimension<3;dimension++) {
        /* Scan all the residents & immigrants to list moved-out atoms */
        for (i=0; i<num_atoms_local+new_atoms; i++) {
            node_lower = 2*dimension;
            node_higher = 2*dimension+1;

            /* Register a to-be-copied atom in move_queue[kul|kuh][] */
            if (!atom_moved[i]) { /* Don't scan moved-out atoms */
                // Check if this atom moved
                if (atom_did_change_node(&positions[3*i],node_lower)) move_queue[node_lower][ ++move_queue[node_lower][0] ] = i;
                else if (atom_did_change_node(&positions[3*i],node_higher)) move_queue[node_higher][ ++move_queue[node_higher][0] ] = i;
            }
        }

        /* Message passing with neighbor nodes----------------------------*/
        com1 = MPI_Wtime();

        /* Loop over the lower & higher directions------------------------*/
        for (j=0; j<2; j++) {
            local_node_id=2*dimension+j;
            node_id = neighbor_nodes[local_node_id]; /* Neighbor node ID */
            /* Send atom-number information---------------------------------*/

            num_send = move_queue[local_node_id][0]; /* # of atoms to-be-sent */

            /* Even node: send & recv*/
            /* Odd node: recv & send */
            if (my_parity[dimension] == 0) {
                MPI_Send(&num_send,1,MPI_INT,node_id,110,MPI_COMM_WORLD);
                MPI_Recv(&num_receive,1,MPI_INT,MPI_ANY_SOURCE,110,MPI_COMM_WORLD,&status);
            }
            else if (my_parity[dimension] == 1) {
                MPI_Recv(&num_receive,1,MPI_INT,MPI_ANY_SOURCE,110, MPI_COMM_WORLD,&status);
                MPI_Send(&num_send,1,MPI_INT,node_id,110,MPI_COMM_WORLD);
            }
            /* Single layer: Exchange information with myself */
            else num_receive = num_send;
            /* Send & receive information on boundary atoms-----------------*/

            /* Message buffering */
            for (i=1; i<=num_send; i++) {
                for (a=0; a<3; a++) {
                    /* Shift the coordinate origin */
                    mpi_send_buffer[6*(i-1)    + a] = positions [3*move_queue[local_node_id][i]+a] - shift_vector[local_node_id][a];
                    mpi_send_buffer[6*(i-1)+ 3 + a] = velocities[3*move_queue[local_node_id][i]+a];
                    atom_moved[ move_queue[local_node_id][i] ] = true;
                }
            }

            /* Even node: send & recv, if not empty */
            /* Odd node: recv & send, if not empty */
            if (my_parity[dimension] == 0) {
                MPI_Send(mpi_send_buffer,6*num_send,MPI_DOUBLE,node_id,120,MPI_COMM_WORLD);
                MPI_Recv(mpi_receive_buffer,6*num_receive,MPI_DOUBLE,MPI_ANY_SOURCE,120,MPI_COMM_WORLD,&status);
            }
            else if (my_parity[dimension] == 1) {
                MPI_Recv(mpi_receive_buffer,6*num_receive,MPI_DOUBLE,MPI_ANY_SOURCE,120,MPI_COMM_WORLD,&status);
                MPI_Send(mpi_send_buffer,6*num_send,MPI_DOUBLE,node_id,120,MPI_COMM_WORLD);
            }
            /* Single layer: Exchange information with myself */
            else for (i=0; i<6*num_receive; i++) mpi_receive_buffer[i] = mpi_send_buffer[i];

            /* Message storing */
            for (i=0; i<num_receive; i++) {
                for (a=0; a<3; a++) {
                    positions [3*(num_atoms_local+new_atoms+i) + a] = mpi_receive_buffer[6*i   + a];
                    velocities[3*(num_atoms_local+new_atoms+i) + a] = mpi_receive_buffer[6*i+3 + a];
                    atom_moved[num_atoms_local+new_atoms+i] = false;
                }
            }

            /* Increment the # of new atoms */
            new_atoms += num_receive;

            /* Internode synchronization */
            MPI_Barrier(MPI_COMM_WORLD);
        } /* Endfor lower & higher directions, j */

        // comt+=MPI_Wtime()-com1;
    }
    int ipt = 0;
    for (i=0; i<num_atoms_local+new_atoms; i++) {
        if (!atom_moved[i]) {
            for (a=0; a<3; a++) {
                positions [3*ipt+a] = positions [3*i+a];
                velocities[3*ipt+a] = velocities[3*i+a];
            }

            ipt++;
        }
    }

    /* Update the compressed # of resident atoms */
    num_atoms_local = ipt;

}

void System::half_kick() {

}

void System::full_kick() {

}

void System::calculate_accelerations() {

}

void System::move() {
    for(n=0;n<num_atoms_local;n++) {
        for(a=0;a<3;a++) {
            positions[3*n+a] += velocities[3*n+a]*dt;
        }

        atom_moved[n] = false;
    }
}

void System::step() {
    if(myid==0 && !(steps % 100) ) cout << steps << endl;
    move();
    mpi_move();
    mpi_copy();
    steps++;
}
