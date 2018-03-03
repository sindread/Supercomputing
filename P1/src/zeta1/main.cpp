#include "zeta1.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <limits>

#define M_PI acos(-1.0)

using namespace std;


bool unit_test()
{
	const auto answer_should_be = &zeta1_expected_value_after_3_iterations;

	return *answer_should_be == run_calculation(2, 3);
}

string verification_test(const int maxk)
{
	ostringstream oss;
	auto n = 0;
	for (auto k = 1; k <= maxk; k++)
	{
		n = pow(2, k);
		oss << "n = 2^" << k << " = " << n << " M_PI - pi_n = " << abs(M_PI - run_calculation(2, n)) << endl;
	}

	return oss.str();
}

int main(int argc, char* argv[])
{
	auto n = 3;
	auto argument_number = 1;

	cout.precision(numeric_limits<double>::digits10 + 2);

	if (argc == argument_number)
	{
		cout << fixed << "Calculate pi for zeta1 function: " << run_calculation(2, n) << endl;
	}

	else
	{
		string arg = argv[argument_number];
		if (arg == "-u")
		{
			cout << fixed << "zeta1 unit test result, with n = 3: " << boolalpha << unit_test() << endl;
		}
		else if (arg== "-v")
		{
			cout << fixed << "Running zeta1 verification tests: " << endl;
			cout << verification_test(24);
		}
		else{
			int argument = stoi(arg);
			if((argument & (argument-1)) != 0 && argument != 0){
				cout << "Number of processes need to be power of two";
				return -1;
			}
			cout << fixed << "Running zeta1 with " << arg << " processes." << endl;
			cout << fixed << run_calculation(2, argument) << endl;
		}
	}
	return 0;
}