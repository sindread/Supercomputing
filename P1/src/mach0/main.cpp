#include "mach0.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <limits>
#include <fstream>

#define M_PI acos(-1.0)

using namespace std;

double calculate_pi(const int &n)
{
	return mach0_calculate_pi(n);
}

bool unit_test()
{
	const auto answer_should_be = &mach0_expected_value_after_3_iterations;

	return *answer_should_be == calculate_pi(3);
}

void verification_test(const unsigned int maxk)
{
	ofstream write_to_file("../verification/mach0 verification.txt");
	auto n = 0;

	write_to_file << "mach0 verification tests " << endl;
	for (auto k = 1; k <= maxk; k++)
	{
		n = pow(2, k);
		write_to_file << "n =" << n << ", error: M_PI - pi_" << n << " = " << abs(M_PI - calculate_pi(n)) << endl;
	}
}

int main(int argc, char* argv[])
{
	auto n = 3;
	auto argument_number = 1;

	cout.precision(numeric_limits<double>::digits10 + 2);

	if (argc == argument_number)
	{
		cout << fixed << "Calculate pi for mach0 function: " << endl << calculate_pi(n) << endl;
	}

	else
	{
		string arg = argv[argument_number];
		if (arg == "-u")
		{
			cout << fixed << "mach0 unit test result, with n = 3: " << boolalpha << unit_test() << endl;
		}
		else if (arg== "-v")
		{
			cout << fixed << "Running mach0 verification tests: " << endl;
			verification_test(24);
		}
	}

	return 0;
}