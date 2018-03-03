#include "mach1.h"
#include <cmath>
#include <iostream>
#include <mpi.h>

using namespace std;


double arctan(const int& nOP, const int& n, const double& x, const int& rank)
{
	auto s = 0.0;

	for (auto i = rank + 1; i <= n; i += nOP)
	{
		const auto part1 = pow(-1, i - 1);
		const auto part3 = (2 * i) - 1;
		const auto part2 = pow(x, part3);
		s += part1 * (part2 / part3);
	}

	return s;
}


void mach1_calculate_vi(const int& nOP, const int& n, const int& rank, double& answer)
{
	const auto x1 = double(1)/5;
	const auto x2 = double(1)/239;

	// pi = 4 * (4 * arctan(1/5) - arctan (1/239))
	answer = (4 * arctan(nOP, n, x1, rank) - arctan(nOP, n, x2, rank));
}



double mach1_calculate_pi(int* argc, char*** argv, const int& n) 
{
	int nOP, rank;
	// int lengthForRank[nOP];
	// int length = n / nOP;
	// int mod = n % nOP;
	
	// for(int i = 0; i < nOP; i++){
	// 	int temp = length;
	// 	if(i < mod){
	// 		temp += 1;
	// 	}
	// 	cout << "Length for rank " << i << temp << endl;
	// 	lengthForRank[i] = temp;
	// }
	
	MPI_Init(argc,argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nOP);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	double startTime;
	if(rank == 0) startTime = MPI_Wtime();

	//double answers[lengthForRank[rank]];
	double answer = 0.0;

	mach1_calculate_vi(nOP, n, rank, answer);
	cout << nOP << " " << n << " " << rank << endl;
	cout << answer << " for rank " << rank << endl;
	//double myPI = 0.0;

	//sumVector(answers, lengthForRank[rank], myPI);
	
	double pi;
	MPI_Reduce(&answer, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(rank == 0){
		cout << "PI:" << pi << endl; 
		pi = 4*pi;
	}

	MPI_Finalize();
	
	return pi;
}
