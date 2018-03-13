#include "zeta0.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>

#define M_PI acos(-1.0)

using namespace std;

double calculate_pi(const int &n)
{
	return zeta0_calculate_pi(n);
}

bool unit_test()
{
	const auto answer_should_be = &zeta0_expected_value_after_3_iterations;

	return *answer_should_be == calculate_pi(3);
}

void verification_test(const int maxk)
{
	ofstream write_to_file("../verification/zeta0 verification.txt");
	auto n = 0;

	write_to_file << "zeta0 verification tests " << endl;
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
		cout << fixed << "Calculate pi for zeta0 function: " << calculate_pi(n) << endl;
	}

	else
	{
		string arg = argv[argument_number];
		if (arg == "-u")
		{
			cout << fixed << "zeta0 unit test result, with n = 3: " << boolalpha << unit_test() << endl;
		}
		else if (arg== "-v")
		{
			cout << fixed << "Running zeta0 verification tests, writing to file" << endl;
			verification_test(24);
		}
	}

	return 0;
}