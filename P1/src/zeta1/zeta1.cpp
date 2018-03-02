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

void sumVector(const double* vector, const int& length, double& sum){
	for(int i = 0; i < length; i++){
		sum += vector[i];	
	}
}


double run_calculation(const int &numberOfProcessors, const int &numberOfIntervals){
	int nOP, rank;
	nOP = numberOfProcessors;
	int lengthForRank[nOP];
	int length = numberOfIntervals / nOP;
	int mod = numberOfIntervals % nOP;
	
	for(int i = 0; i < nOP; i++){
		int temp = length;
		if(i < mod){
			temp += 1;
		}
		cout << "Length for rank " << i << temp << endl;
		lengthForRank[i] = temp;
	}
	
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &nOP);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	double startTime;
	if(rank == 0) startTime = MPI_Wtime();

	double answers[lengthForRank[rank]];

	zeta1_calculate_vi(nOP, rank, numberOfIntervals, answers);

	double myPI = 0.0;

	sumVector(answers, lengthForRank[rank], myPI);
	
	double pi;
	MPI_Reduce(&myPI, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(rank == 0){
		cout << "PI:" << pi << endl; 
		pi = sqrt(pi * 6);
	}

	MPI_Finalize();
	
	return pi;
	
}



