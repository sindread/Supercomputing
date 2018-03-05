#include "zeta2.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <limits>

#define M_PI acos(-1.0)

using namespace std;

double calculate_pi(const int &n)
{
	return zeta2_calculate_pi(n);
}

bool unit_test()
{
	const auto answer_should_be = &zeta2_expected_value_after_3_iterations;

	return *answer_should_be == calculate_pi(3);
}

string verification_test(const int maxk)
{
	ostringstream oss;
	auto n = 0;
	for (auto k = 1; k <= maxk; k++)
	{
		n = pow(2, k);

		auto answer = abs(M_PI - calculate_pi(n));

		oss << "n =" << n << ", error: M_PI - pi_n = " << answer << endl;
	}

	return oss.str();
}

int main(int argc, char* argv[])
{
	auto n = 10;
	auto argument_number = 1;

	cout.precision(numeric_limits<double>::digits10 + 2);
	
	if (argc == argument_number) {
		
		auto answer = calculate_pi(n);
		cout << fixed << "Calculate pi for zeta2 function: " << answer << endl;
	
	} else {
		string arg = argv[argument_number];
		if (arg == "-u")
		{
			auto boolalpha = unit_test();
			cout << fixed << "zeta2 unit test result, with n = 3: " << boolalpha << endl;
		}
		else {
			// auto argument = stoi(arg);
			// if ( (argument & (argument-1) ) != 0 && argument != 0) {
			// 	cout << fixed << "Number of processes need to be power of two";
			// 	return -1;
			// }

			auto answer = calculate_pi(n);

			cout << fixed << "Running zeta2 with " << arg << " processes." << endl;
			cout << fixed << answer << endl;
		}
	}


	return 0;
}