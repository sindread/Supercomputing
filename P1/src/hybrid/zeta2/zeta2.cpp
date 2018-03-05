#include "zeta2.h"
#include <cmath>
#include <omp.h>
#include <chrono>
#include <iostream>
#include <mpi.h>

using namespace std;

void zeta2_calculate_vi(const int &numberOfProcesses, const int &rank, const int &n, double &sumPiParts)
{
	#pragma omp parallel for reduction (+:sumPiParts)
	for (int i = rank + 1; i <= n; i += numberOfProcesses)
	{
		cout << omp_get_thread_num() << endl;
		sumPiParts += 1.0 / (i*i);
	}
}

void sumVector(const double* vector, const int& length, double& sum){
	for(int i = 0; i < length; i++){
		sum += vector[i];	
	}
}

double zeta2_calculate_pi(const int &rank, const int &numberOfProcessors, const int &numberOfIntervals){
	int length = numberOfIntervals / numberOfProcessors + (numberOfIntervals % numberOfProcessors == 0 ? 0 : 1);

	double startTime;
	if(rank == 0){
		startTime = MPI_Wtime();
	} 

	// double piParts[length];
	// double sumPiParts = 0.0; 
	
	// zeta2_calculate_vi(numberOfProcessors, rank, numberOfIntervals, piParts);
	// sumVector(piParts, length, sumPiParts);

	double sumPiParts = 0.0; 
	zeta2_calculate_vi(numberOfProcessors, rank, numberOfIntervals, sumPiParts);
	
	double pi; 
	MPI_Reduce(&sumPiParts, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if(rank == 0){
		pi = sqrt(pi * 6);

		double endTime = MPI_Wtime();

		cout << "Main process used " << (endTime - startTime) << endl;
	} else {
		cout << "Helper process " << rank << " done."<< endl;  
	}
	
	return pi;
}



