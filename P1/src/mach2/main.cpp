#include "mach2.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <limits>
#include <mpi.h>

#define M_PI acos(-1.0)

using namespace std;

double calculate_pi(const int &rank, const int &numberOfprocesses, const int &n)
{
	return mach2_calculate_pi(rank , numberOfprocesses, n);
}

bool unit_test(int &rank)
{
	const auto answer_should_be = &mach2_expected_value_after_3_iterations;

	return *answer_should_be == calculate_pi(rank, 2, 3);
}

string verification_test(int &rank, const int maxk)
{
	ostringstream oss;
	auto n = 0;
	for (auto k = 1; k <= maxk; k++)
	{
		n = pow(2, k);

		auto answer = abs(M_PI - calculate_pi(rank, 2, n));

		if (rank == 0){
			oss << "n =" << n << ", error: M_PI - pi_n = " << answer << endl;
		}
	}

	return oss.str();
}

int main(int argc, char* argv[])
{
	auto n = 100;
	auto argument_number = 1;
	auto numberOfprocesses = 2;
	int rank;

	MPI_Init(&argc , &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfprocesses);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		cout.precision(numeric_limits<double>::digits10 + 2);
	}
	
	if (argc == argument_number) {
		// Parallel calculation
		auto answer = calculate_pi(rank, 2, n);

		// Seriel print:
		if (rank == 0){
			cout << fixed << "Calculate pi for mach2 function: " << answer << endl;
		}
	} else {
		string arg = argv[argument_number];
		if (arg == "-u")
		{
			// Parallel calculation
			auto boolalpha = unit_test(rank);

			// Seriel print:
			if (rank == 0){
				cout << fixed << "mach2 unit test result, with n = 3: " << boolalpha << endl;
			}
		}
		else {
			auto argument = stoi(arg);
			// Seriel start:
			if (rank == 0){
				if ( (argument & (argument-1) ) != 0 && argument != 0) {
					cout << fixed << "Number of processes need to be power of two";
					return -1;
				}
			}

			// Parallel calculation
			auto answer = calculate_pi(rank, 2, argument);

			// Seriel print:
			if (rank == 0){
				cout << fixed << "Running mach2 with " << arg << " processes." << endl;
				cout << fixed << answer << endl;
			}
		}
	}

	MPI_Finalize();
	return 0;
}