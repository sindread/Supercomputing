#include "zeta1.h"
#include <iostream>
#include <mpi.h>
#include <cmath>

using namespace std;

int main(int argc, char* argv[])
{
	int numberOfProcesses, rank, p = 10;

	MPI_Init(NULL , NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if ( (numberOfProcesses & (numberOfProcesses-1) ) != 0 && numberOfProcesses != 0) {
		if (rank = 0){
			cout << fixed << "Number of processes need to be power of two"<< endl;
		}
		
		MPI_Finalize();
		return -1;
	}

	if (argc > 1){
		string arg = argv[1];
		if (arg =="-v"){
			p = 24;
		}
	} else {
		p = -1;
	}

	if (rank == 0){
		int n = 3;
		master_task(n, numberOfProcesses);
			
	} else {
		slave_task(rank, numberOfProcesses);
	}
	
	MPI_Finalize();
	return 0;
}