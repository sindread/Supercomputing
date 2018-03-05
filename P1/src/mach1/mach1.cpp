#include "mach1.h"
#include <cmath>
#include <iostream>
#include <mpi.h>

using namespace std;

void master_init(int argc, char* argv[], int &n){
	if (argc > 1){
		string arg = argv[1];
		if (arg =="-n"){
			n = stoi(argv[2]);
		}
	} else {
		n = 3;
	}
}

void master_task(const int &n, const int &numberOfProcesses){
	double vi5[n], vi239[n];
	double start, end;

	start = MPI_Wtime();

	vi_parts(n, vi5, vi239);

	int lengthForRank[numberOfProcesses];
	length_of_work(lengthForRank, n, numberOfProcesses);
	
	MPI_Bcast(&lengthForRank, numberOfProcesses-1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&vi5, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&vi239, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	double piSum, pi;
	MPI_Reduce(&piSum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	end = MPI_Wtime();

	cout << "Pi is with mach1, with " << n << " iterations: Pi = " << pi <<  endl;
	cout << "Error(PI-pi_" << n << "): E  = " << M_PI-pi <<  endl;
	cout << "Runtime: Time = " <<  (end-start)*1000 << "ms" << endl;
}

void slave_task(int &rank, int &numberOfProcesses){
	int n;
	int slaves = numberOfProcesses-1;
	int workeRank = rank -1;
	int lengthForRank[workeRank];
	
	MPI_Bcast(&lengthForRank, slaves, MPI_INT, 0, MPI_COMM_WORLD);
	sumVector(lengthForRank, slaves, n);

	double vi5[n], vi239[n];
	MPI_Bcast(&vi5, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&vi239, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int length = lengthForRank[workeRank];
	int index = 0;
	sumVector(lengthForRank, workeRank, index);

	double piPart5, piPart239;
	sumVector(vi5, length, piPart5);
	sumVector(vi239, length, piPart239);

	// pi = 4 * (4 * arctan(1/5) - arctan (1/239))
	double piSum = 4 * (4 * piPart5 - piPart239);

	MPI_Reduce(&piSum, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
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

	return part1 * (part2 / part3);
}

void sumVector(const double* vector, const int& length, double& sum){
	for(int i = 0; i < length; i++){
		sum += vector[i];	
	}
}

void sumVector(const int* vector, const int& length, int& sum){
	for(int i = 0; i < length; i++){
		sum += vector[i];	
	}
}