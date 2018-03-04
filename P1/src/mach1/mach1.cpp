#include "mach1.h"
#include <cmath>
#include <iostream>
#include <mpi.h>

using namespace std;

void master_init(int argc, char* argv[], int &n){
	if (argc > 1){
		string arg = argv[1];
		else if (arg =="-n"){
			n = stoi(argv[2]);
		}
	} else {
		n = 3;
	}
}

void master_task(const int &n, const int &numberOfProcesses){
	double vi5[n], vi239[n];
	double piSum, piSumPpart, pi, start, end;

	start = MPI_Wtime();

	vi_parts(n, vi5, vi239);

	int lengthForRank[numberOfProcesses];
	length_of_work(lengthForRank, n, numberOfProcesses)

	int index = 0;
	for (auto i = 1; i < numberOfProcesses; i++){
		MPI_Send(&lengthForRank[i-1], 1, MPI_INT, i, TAG_LENGTH, MPI_COMM_WORLD);
		MPI_Send(&vi5[index], lengthForRank[i-1], MPI_DOUBLE, i, TAG_VPARTS_5, MPI_COMM_WORLD);
		MPI_Send(&vi239[index], lengthForRank[i-1], MPI_DOUBLE, i, TAG_VPARTS_239, MPI_COMM_WORLD);
		index += lengthForRank[i-1];
	}
	
	int sources = 1;
	while (sources < numberOfProcesses){
		MPI_Recv(&piSumPpart, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		piSum += piSumPpart;
		sources++;
	}

	pi = piSum;

	end = MPI_Wtime();

	cout << "Pi is with mach1, with " << n << " iterations: Pi = " << pi <<  endl;
	cout << "Error(PI-pi_" << n << "): E  = " << M_PI-pi <<  endl;
	cout << "Runtime: Time =" <<  (end-start)*1000 << "ms" << endl;
}

void slave_task(int &rank, int &numberOfProcesses){
	int length;
	MPI_Recv(&length, 1, MPI_INT, 0, TAG_LENGTH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	double vi_parts5[length];
	double vi_parts239[length];
	MPI_Recv(&vi_parts5, length, MPI_DOUBLE, 0, TAG_VPARTS_5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(&vi_parts239, length, MPI_DOUBLE, 0, TAG_VPARTS_239, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	double piPart5, piPart239;
	sumVector(vi_parts5, length, piPart5);
	sumVector(vi_parts239, length, piPart239);

	// pi = 4 * (4 * arctan(1/5) - arctan (1/239))

	double piPart = 4 * (4 * piPart5 - piPart239);

	MPI_Send(&piPart, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
}

void vi_parts(const int& n, double* vi5, double* vi239)
{
	const auto x1 = double(1)/5;
	const auto x2 = double(1)/239;

	for (int i = 1; i <= n; i++){
		vi5[i-1] = arctan_part(i, x1);
		vi239[i-1] = arctan_part(i, x2);
	}
}

double arctan_part(const int& i, const double& x)
{
	const auto part1 = pow(-1, i - 1);
	const auto part3 = (2 * i) - 1;
	const auto part2 = pow(x, part3); 

	return part1 * (part2 / part3);;
}

void sumVector(const double* vector, const int& length, double& sum){
	for(int i = 0; i < length; i++){
		sum += vector[i];	
	}
}
