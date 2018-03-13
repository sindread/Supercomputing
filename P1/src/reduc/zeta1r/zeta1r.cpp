#include "zeta1r.h"
#include <cmath>
#include <mpi.h>
#include <stdio.h>

#define M_PI acos(-1.0)

using namespace std;

void master_task(const int &n, const int &numberOfProcesses){
	double vi[n];
	double start, end;

	start = MPI_Wtime();

	vi_parts(n, vi);

	int lengthForRank[numberOfProcesses];
	length_of_work(lengthForRank, n, numberOfProcesses);

	MPI_Bcast(&lengthForRank, numberOfProcesses-1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&vi, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	//Reduction and broadcast with MPI
	double sum, pi;
	MPI_Allreduce(&sum, &pi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	pi = sqrt(6*pi);
	
	//Manual global reduction
	double sum2;
	globalReduce(numberOfProcesses, sum2, lengthForRank);

	MPI_Bcast(&sum2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	end = MPI_Wtime();
	
	printf("Pi is with zeta1, with ", n," iterations: Pi = ", pi );
	printf("Error(PI-pi_" ,n ,"): E  = " ,M_PI-pi );
	printf("Runtime: Time = " , (end-start)*1000 ,"ms" );
}

void slave_task(int &rank, int &numberOfProcesses){
	int n;
	int slaves = numberOfProcesses-1;
	int workerRank = rank -1;
	int lengthForRank[workerRank];

	MPI_Bcast(&lengthForRank, slaves, MPI_INT, 0, MPI_COMM_WORLD);
	sumVector(lengthForRank, slaves, n);

	double vi[n];
	MPI_Bcast(&vi, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int length = lengthForRank[workerRank];
	int index = 0;
	sumVector(lengthForRank, workerRank, index);

	double partSum = 0.0;
	sumVector(&vi[index], length, partSum);

	//Reduction and broadcast with MPI
	double sum;
	MPI_Allreduce(&partSum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	sum = sqrt(6*sum);
	printf("Sum in slave number " ,rank ," after global reduction (MPI): " ,sum );
	
	// //Manual global reduction
	MPI_Send(&vi[index], length, MPI_DOUBLE, 0, TAG_PARTSUM, MPI_COMM_WORLD);
	double sum2 = 0.0;

	MPI_Bcast(&sum2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	printf("Sum in slave number " ,rank ," after global reduction (manual): " ,sum2 );
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

void globalReduce(const int &numberOfProcesses, double &sum, const int *lengthForProcesses){
	//MPI_Allgather(&myArray, length, MPI_DOUBLE, array&, )
	for (int i = 1; i < numberOfProcesses; i++){
		int length = lengthForProcesses[i-1];
		double array[length];
		MPI_Recv(&array, length, MPI_DOUBLE, i, TAG_PARTSUM, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		sumVector(array, length, sum);
	}
	sum = sqrt(sum*6);
}