#include "mach0.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <limits>

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

string verification_test(const unsigned int maxk)
{
	ostringstream oss;
	auto n = 0;
	for (auto k = 1; k <= maxk; k++)
	{
		n = pow(2, k);
		oss << "n = 2^" << k << " = " << n << " M_PI - pi_n = " << abs(M_PI - calculate_pi(n)) << endl;
	}

	return oss.str();
}

int main(int argc, char* argv[])
{
	auto n = 3;
	cout.precision(numeric_limits<double>::digits10 + 2);

	if (argc == 0)
	{
		cout << fixed << "Calculate pi for zeta0 function: " << calculate_pi(n) << endl;
	}

	else
	{
		if (argv[0] == "-u")
		{
			cout << fixed << "Zeta0 unit test result, with n = 3: " << boolalpha << unit_test() << endl;
		}
		else if (argv[0] == "-v")
		{
			cout << fixed << "Running verification tests: " << endl;
			cout << verification_test(24);
		}
	}

	cin.get();
	return 0;
}