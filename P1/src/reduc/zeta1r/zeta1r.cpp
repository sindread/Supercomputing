#include "zeta1r.h"
#include <cmath>
#include <mpi.h>

#define M_PI acos(-1.0)

using namespace std;

double sum;

void master_init(int argc, char* argv[], int &n){
	if (argc > 1){
		string arg = argv[1];
		if (arg =="-n"){
			n = stoi(argv[2]);
		}
	} else {
		n = 100000;
	}
}

void master_task(const int &n, const int &numberOfProcesses){
	double vi[n];
	double start, end;

	start = MPI_Wtime();

	vi_parts(n, vi);

	int lengthForRank[numberOfProcesses];
	length_of_work(lengthForRank, n, numberOfProcesses);

	int index = 0;
	for (auto i = 1; i < numberOfProcesses; i++){
		MPI_Send(&lengthForRank[i-1], 1, MPI_INT, i, TAG_LENGTH, MPI_COMM_WORLD);
		MPI_Send(&vi[index], lengthForRank[i-1], MPI_DOUBLE, i, TAG_VPARTS, MPI_COMM_WORLD);
		index += lengthForRank[i-1];
	}
	
	//Reduction and broadcast with MPI
	double pi;
	MPI_Allreduce(&sum, &pi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	pi = sqrt(6*pi);
	
	//Manual global reduction
	double sum2;
	globalReduce(numberOfProcesses, sum2, lengthForRank);
	MPI_Bcast(&sum2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	cout << sum2 << endl;

	end = MPI_Wtime();
	
	cout << "Pi is with zeta1, with " << n << " iterations: Pi = " << pi <<  endl;
	cout << "Error(PI-pi_" << n << "): E  = " << M_PI-pi <<  endl;
	cout << "Runtime: Time = " <<  (end-start)*1000 << "ms" << endl;
}

void slave_task(int &rank, int &numberOfProcesses){
	int length;
	MPI_Recv(&length, 1, MPI_INT, 0, TAG_LENGTH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	double vi_parts[length];
	MPI_Recv(&vi_parts, length, MPI_DOUBLE, 0, TAG_VPARTS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	double partSum = 0.0;
	sumVector(vi_parts, length, partSum);

	//Reduction and broadcast with MPI
	MPI_Allreduce(&partSum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	cout << "Sum in slave number " << rank << ": " << sum << endl;

	//Manual global reduction
	MPI_Send(&vi_parts, length, MPI_DOUBLE, 0, TAG_PARTSUM, MPI_COMM_WORLD);
	double sum2 = 0.0;
	MPI_Bcast(&sum2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	cout << "Sum from rank " << rank << ": " << sum2 << endl;
	

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