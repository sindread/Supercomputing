#include "zeta1.h"
#include <cmath>
#include <mpi.h>

#define M_PI acos(-1.0)

using namespace std;

void master_init(int argc, char* argv[], int &n){
	if (argc > 1){
		string arg = argv[1];
		if (arg =="-n"){
			n = stoi(argv[2]);
		}
	} else {
		n = 100;
	}
}

void master_task(const int &n, const int &numberOfProcesses){
	double vi[n];
	double start, end;

	start = MPI_Wtime();

	vi_parts(n, vi);

	int lengthForRank[numberOfProcesses];
	length_of_work(lengthForRank, n, numberOfProcesses)

	int index = 0;
	for (auto i = 1; i < numberOfProcesses; i++){
		MPI_Send(&lengthForRank[i-1], 1, MPI_INT, i, TAG_LENGTH, MPI_COMM_WORLD);
		MPI_Send(&vi[index], lengthForRank[i-1], MPI_DOUBLE, i, TAG_VPARTS, MPI_COMM_WORLD);
		index += lengthForRank[i-1];
	}
	
	int sources = 1;
	double piSum, piSumPpart, pi;
	while (sources < numberOfProcesses){
		MPI_Recv(&piSumPpart, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		piSum += piSumPpart;
		sources++;
	}

	pi = sqrt(6*piSum);
	end = MPI_Wtime();
	
	cout << "Pi is with mach1, with " << n << " iterations: Pi = " << pi <<  endl;
	cout << "Error(PI-pi_" << n << "): E  = " << M_PI-pi <<  endl;
	cout << "Runtime: Time =" <<  (end-start)*1000 << "ms" << endl;
}

void slave_task(int &rank, int &numberOfProcesses){
	int length;
	MPI_Recv(&length, 1, MPI_INT, 0, TAG_LENGTH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	double vi_parts[length];
	MPI_Recv(&vi_parts, length, MPI_DOUBLE, 0, TAG_VPARTS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	double piPart;
	sumVector(vi_parts, length, piPart);

	MPI_Send(&piPart, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
}

void vi_parts(const int &n, double* vi)
{
	for (int i = 1; i < n; i++)
	{
		 vi[i-1] = 1.0 / pow(i,2);
	}
}

void length_of_work(int* lengthForRank, const int &n, const int &numberOfProcesses){
	int length = n / (numberOfProcesses-1);
	int rest = n % (numberOfProcesses-1);
	
	for(int i = 0; i < numberOfProcesses-1; i++){
		lengthForRank[i] = length;

		if(i < rest){
			lengthForRank[i]++;
		}
	}
}

void sumVector(const double* vector, const int& length, double& sum){
	for(int i = 0; i < length; i++){
		sum += vector[i];	
	}
}