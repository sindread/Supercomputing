#include "zeta1.h"
#include <cmath>
#include <mpi.h>

using namespace std;

void master_init(int argc, char* argv[], int &n, bool &isUnitTest){
	if (argc == 2){
		n = stoi(argv[1]);
	} else if (argc == 3) {
		auto arg = argv[2];
		if (arg == "-u"){
			isUnitTest = true;
			n = 3;
		}
	} else {
		n = 4;
	}
}

void master_task(const int &n, const int &numberOfProcesses){
	double* vi = vi_parts(n);

	int lengthForRank[numberOfProcesses];
	int length = n / (numberOfProcesses-1);
	int rest = n % (numberOfProcesses-1);
	
	for(int i = 0; i < numberOfProcesses-1; i++){
		lengthForRank[i] = length;

		if(i < rest){
			lengthForRank[i]++;
		}
	}

	int counter = 0;

	for (auto i = 1; i < numberOfProcesses; i++){
		if (rest > 0){
			MPI_Send(&lengthForRank[i-1], 1, MPI_INT, i, TAG_LENGTH, MPI_COMM_WORLD);
			MPI_Send(&vi, lengthForRank[i-1], MPI_DOUBLE, i, TAG_VPARTS, MPI_COMM_WORLD);
			vi = vi + lengthForRank[i-1];
		}
	}

	double sum = 0.0;
	double part;
	int sources = 1;
	while (sources < numberOfProcesses){
		MPI_Recv(&part, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		sum += part;
		sources++;
	}

	double pi = sqrt(6*sum);

	cout << pi << endl;
}

void slave_task(int &rank, int &numberOfProcesses){
	int length;
	MPI_Recv(&length, 1, MPI_INT, 0, TAG_LENGTH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	double vi_parts[length];
	MPI_Recv(vi_parts, length, MPI_DOUBLE, 0, TAG_VPARTS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	double piPart;
	sumVector(vi_parts, length, piPart);

	MPI_Send(&piPart, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
}

double* vi_parts(const int &n)
{
	double parts[n];
	for (int i = 1; i < n; i++)
	{
		 parts[i-1] = 1.0 / pow(i,2);
	}

	return parts;
}

void sumVector(const double* vector, const int& length, double& sum){
	for(int i = 0; i < length; i++){
		sum += vector[i];	
	}
}