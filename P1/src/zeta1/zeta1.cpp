#include "zeta1.h"
#include <cmath>
#include <mpi.h>

#define M_PI acos(-1.0)

using namespace std;

void master_task(const int &n, const int &numberOfProcesses){
	double vi[n];
	double start, end;

	start = MPI_Wtime();

	vi_parts(n, vi);

	int lengthForRank[numberOfProcesses-1];
	length_of_work(lengthForRank, n, numberOfProcesses);

	MPI_Bcast(&lengthForRank, numberOfProcesses-1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Bcast(&vi, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	double piSum = 0, pi = 0;
	MPI_Reduce(&piSum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	pi = sqrt(6*pi);
	end = MPI_Wtime();
	
	cout << "Pi is with zeta1, with " << n << " iterations: Pi = " << pi <<  endl;
	cout << "Error(PI-pi_" << n << "): E  = " << abs(M_PI-pi) <<  endl;
	cout << "Runtime: Time = " <<  (end-start)*1000 << "ms" << endl;
}

void slave_task(int &rank, int &numberOfProcesses){
	int n;
	int slaves = numberOfProcesses-1;
	int WorkerRank = rank -1;
	int lengthForRank[WorkerRank];
	
	MPI_Bcast(lengthForRank, slaves, MPI_INT, 0, MPI_COMM_WORLD);
	sumVector(lengthForRank, slaves, n);

	double vi[n];
	MPI_Bcast(vi, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	int length = lengthForRank[WorkerRank];
	int index = 0;
	sumVector(lengthForRank, WorkerRank, index);

	double piSum=0;
	sumVector(&vi[index], length ,piSum);
	
	MPI_Reduce(&piSum, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

void vi_parts(const int &n, double* vi)
{
	for (int i = 1; i <= n; i++)
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

void sumVector(const int* vector, const int& length, int& sum){
	for(int i = 0; i < length; i++){
		sum += vector[i];	
	}
}