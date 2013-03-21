#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#define MAX_ATOM_NUM 100000
char atom_type_string[][5] = {"Ar ", "H "};

using namespace std;

int main(int args, char *argv[]) {
	if(args < 2) {
		cout << "Please specify the number of cpus and timesteps." << endl;
		return 0;
	}
	int cpus = atoi(argv[1]);
	double *positions = new double[9*MAX_ATOM_NUM];
	unsigned long *atom_type = new unsigned long[MAX_ATOM_NUM];

	ofstream file ("state.xyz", ios::out);
	
	ifstream **state_files = new ifstream*[cpus];
	for(int cpu=0;cpu<cpus;cpu++) {
		char *filename = new char[100];
		sprintf(filename,"release/state_files/state%04d.bin",cpu);
		state_files[cpu] = new ifstream(filename,ios::in | ios::binary);
	}
	cout << cpus << " state files opened." << endl;
	
	int num_particles = 0;
	for(int cpu=0;cpu<cpus;cpu++) { 
		int N;
		state_files[cpu]->read(reinterpret_cast<char*>(&N),sizeof(int));
		state_files[cpu]->read(reinterpret_cast<char*>(&positions[6*num_particles]),6*N*sizeof(double));
		state_files[cpu]->read(reinterpret_cast<char*>(&atom_type[num_particles]),N*sizeof(unsigned long));
		num_particles += N;
	}

	file << num_particles << endl;
	file << "sup" << endl;
	for(int n=0;n<num_particles;n++) {
    	file << atom_type_string[atom_type[n]] << positions[6*n+0] << " " << positions[6*n+1] << " " << positions[6*n+2] << endl;
    }

	file.close();
	for(int cpu=0;cpu<cpus;cpu++) {
		state_files[cpu]->close();
	}
	cout << "XYZ-file created." << endl;

	return 0;
}