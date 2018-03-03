#include "mach1.h"
#include <cmath>
#include <iostream>
#include <mpi.h>

using namespace std;


double arctan(const int& numberOfProcessors, const int& n, const double& x, const int& rank)
{
	auto s = 0.0;

	for (auto i = rank + 1; i <= n; i += numberOfProcessors)
	{
		const auto part1 = pow(-1, i - 1);
		const auto part3 = (2 * i) - 1;
		const auto part2 = pow(x, part3);
		s += part1 * (part2 / part3);
	}

	return s;
}

void mach1_calculate_vi(const int& numberOfProcessors, const int& n, const int& rank, double& answer)
{
	const auto x1 = double(1)/5;
	const auto x2 = double(1)/239;

	// pi = 4 * (4 * arctan(1/5) - arctan (1/239))
	answer = (4 * arctan(numberOfProcessors, n, x1, rank) - arctan(numberOfProcessors, n, x2, rank));
}

double mach1_calculate_pi(const int &rank, const int &numberOfProcessors, const int &numberOfIntervals)
{
	double startTime;
	if(rank == 0) startTime = MPI_Wtime();

	double answer = 0.0;
	mach1_calculate_vi(numberOfProcessors, numberOfIntervals, rank, answer);
	
	double pi;
	MPI_Reduce(&answer, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if(rank == 0){
		pi = 4*pi;

		double endTime = MPI_Wtime();

		cout << "Main prosess used " << (endTime - startTime) << endl;
	} else {
		cout << "Helper prosess " << rank << " done."<< endl;  
	}
	
	return pi;
}
