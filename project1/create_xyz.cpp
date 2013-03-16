#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

int main(int args, char *argv[]) {
	if(args < 2) {
		cout << "Please specify the number of cpus and timesteps." << endl;
		return 0;
	}
	int cpus = atoi(argv[1]);
	double *positions = new double[9*1000000];
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
		state_files[cpu]->read(reinterpret_cast<char*>(&positions[9*num_particles]),9*N*sizeof(double));
		num_particles += N;
	}

	file << num_particles << endl;
	file << "sup" << endl;
	for(int n=0;n<num_particles;n++) {
    	file << "H " << positions[9*n+0] << " " << positions[9*n+1] << " " << positions[9*n+2] << endl;
    }

	file.close();
	for(int cpu=0;cpu<cpus;cpu++) {
		state_files[cpu]->close();
	}
	cout << "XYZ-file created." << endl;

	return 0;
}