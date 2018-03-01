#include "zeta1.h"
#include <cmath>
#include <mpi.h>

void zeta1_calculate_vi(const int &numberOfProcesses, const int &rank, const int &n, double* answers)
{
	int j = 0;
	for (int i = rank; i <= n; i += numberOfProcesses)
	{
		const double v_i = 1.0 / (i*i);	
		answers[j] = v_i;
		j++;
	}
	
}


double run_calculation(const int &numberOfIntervals){
	int nOP, rank;
	
	MPI_Init(nullptr, nullptr);
	MPI_Comm_size(MPI_COMM_WORLD, &nOP);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	double startTime;
	if(rank == 0) startTime = MPI_Wtime();
	double answers[numberOfIntervals/nOP];
	zeta1_calculate_vi(nOP, rank, numberOfIntervals, answers);
	
	double pi;
	std::cout << rank << pi << std::endl;
	MPI_Reduce(&answers, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	std::cout << rank << pi << std::endl;
	if(rank == 0){
		pi = sqrt(pi * 6);
	}

	MPI_Finalize();
	
	return pi;
	
}

void sumVector(const double* vector, const int& length, double& sum){
	for(int i = 0; i < length; i++){
		
	}
}


