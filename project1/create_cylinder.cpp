#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#define MAX_ATOM_NUM 100000
using namespace std;

int main(int args, char *argv[]) {
	if(args < 4) {
		cout << "Please specify the number of cpus, radius and cylinder dimension" << endl;
		cout << "./create_cylinder int double [0|1|2]" << endl;
		return 0;
	}

	int cpus = atoi(argv[1]);
	double radius = atof(argv[2]);
	int dimension = atoi(argv[3]);

	double data[6*MAX_ATOM_NUM];
	bool  is_frozen[MAX_ATOM_NUM];
	char *filename = new char[100];
	int num_particles;
	double pos_coeff_squared = 3.405*3.405;
	double radius_squared = radius*radius;
	double r2;

	for(int cpu=0;cpu<cpus;cpu++) { 
		
		sprintf(filename,"release/state_files/state%04d.bin",cpu);
		ifstream state_file(filename,ios::in | ios::binary);

		state_file.read(reinterpret_cast<char*>(&num_particles),sizeof(int));
		state_file.read(reinterpret_cast<char*>(&data),6*num_particles*sizeof(double));
		state_file.read(reinterpret_cast<char*>(&is_frozen),num_particles*sizeof(bool));
		r2 = 0;
		for(int i=0;i<num_particles;i++) {
			for(int a=0;a<3;a++) {
				if(a != dimension) r2 += data[3*i+a]*data[3*i+a]*;
			}
			if(r2 < radius_squared) is_frozen[i] = true;
			else is_frozen[i] = false;
		}
		state_file.close();
		ofstream save_state_file(filename,ios::out | ios::binary);
		save_state_file.write(reinterpret_cast<char*>(&num_particles),sizeof(int));
		save_state_file.write(reinterpret_cast<char*>(&data),6*num_particles*sizeof(double));
		save_state_file.write(reinterpret_cast<char*>(&is_frozen),num_particles*sizeof(bool));
		save_state_file.close();
	}

	cout << "Cylinder created" << endl;
	
	return 0;
}