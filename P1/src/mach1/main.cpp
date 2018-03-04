#include "mach1.h"
#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char* argv[])
{
	int numberOfProcesses, rank;

	MPI_Init(NULL , NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0){
		if ( (numberOfProcesses & (numberOfProcesses-1) ) != 0 && numberOfProcesses != 0) {
			cout << fixed << "Number of processes need to be power of two";
			return -1;
		}

		int n;
		bool isUnitTest = false;
		master_init(argc, argv, n, isUnitTest);
		master_task(n, isUnitTest, numberOfProcesses);
	} else {
		slave_task(rank, numberOfProcesses);
	}

	MPI_Finalize();

	return 0;
}