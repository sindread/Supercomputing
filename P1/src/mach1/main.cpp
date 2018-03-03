#include "mach1.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <limits>

#define M_PI acos(-1.0)

using namespace std;


bool unit_test(int* argc, char** argv[])
{
	const auto answer_should_be = &mach1_expected_value_after_3_iterations;
	char array[1];
	array[0] = 1;
	char* pointer = array;
	auto pointer2 = *(pointer+0);
	cout << pointer << endl;
	cout << pointer2 << endl;
	return false;
	//return *answer_should_be == mach1_calculate_pi(1, pointer2, 3);
}

int main(int argc, char* argv[])
{
	auto n = 3;
	auto argument_number = 1;

	cout.precision(numeric_limits<double>::digits10 + 2);
	
	if (argc == argument_number)
	{
		cout << fixed << "Calculate pi for mach1 function: " << mach1_calculate_pi(&argc, &argv, n) << endl;
	}

	else
	{
		string arg = argv[argument_number];
		if (arg == "-u")
		{
			cout << fixed << "mach1 unit test result, with n = 3: " << boolalpha << unit_test(&argc, &argv) << endl;
		}
		// else if (arg== "-v")
		// {
		// 	cout << fixed << "Running mach1 verification tests: " << endl;
		// 	cout << verification_test(24);
		// }
		else{
			int argument = stoi(arg);
			if((argument & (argument-1)) != 0 && argument != 0){
				cout << "Number of processes need to be power of two";
				return -1;
			}
			cout << fixed << "Running mach1 with " << arg << " processes." << endl;
			cout << fixed << mach1_calculate_pi(&argc, &argv, n) << endl;
		}
	}
	return 0;
}