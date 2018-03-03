#include "zeta1.h"
#include <cmath>
#include <mpi.h>

using namespace std;

void zeta1_calculate_vi(const int &numberOfProcesses, const int &rank, const int &n, double* answers)
{
	int j = 0;
	for (int i = rank + 1; i <= n; i += numberOfProcesses)
	{
		const double v_i = 1.0 / (i*i);	
		answers[j] = v_i;
		j++;
	}
}

void zeta1_calculate_vi(const int &numberOfProcesses, const int &rank, const int &n, double &sumPiParts)
{
	for (int i = rank + 1; i <= n; i += numberOfProcesses)
	{
		sumPiParts += 1.0 / (i*i);
	}
}

void sumVector(const double* vector, const int& length, double& sum){
	for(int i = 0; i < length; i++){
		sum += vector[i];	
	}
}

double zeta1_calculate_pi(int &rank, const int &numberOfProcessors, const int &numberOfIntervals){
	int length = numberOfIntervals / numberOfProcessors + (numberOfIntervals % numberOfProcessors == 0 ? 0 : 1);

	double startTime;
	if(rank == 0){
		startTime = MPI_Wtime();
	} 

	// double piParts[length];
	// double sumPiParts = 0.0; 
	
	// zeta1_calculate_vi(numberOfProcessors, rank, numberOfIntervals, piParts);
	// sumVector(piParts, length, sumPiParts);

	double sumPiParts = 0.0; 
	zeta1_calculate_vi(numberOfProcessors, rank, numberOfIntervals, sumPiParts);
	
	double pi; 
	MPI_Reduce(&sumPiParts, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if(rank == 0){
		pi = sqrt(pi * 6);

		double endTime = MPI_Wtime();

		cout << "Main prosess used " << (endTime - startTime) << endl;
	} else {
		cout << "Helper prosess " << rank << " done."<< endl;  
	}
	
	return pi;
}