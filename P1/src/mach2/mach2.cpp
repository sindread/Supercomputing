#include "mach2.h"
#include <cmath>
#include <iostream>
#include <omp.h>

using namespace std;

double mach2_calculate_pi(const int& n) 
{
	const auto x1 = double(1)/5;
	const auto x2 = double(1)/239;

	// pi = 4 * (4 * arctan(1/5) - arctan (1/239))
	return 4 * (4 * arctan(n, x1) - arctan(n, x2));
}

double arctan(const int& n, const double& x)
{
	auto s = 0.0;

	#pragma omp parallel for reduction (+:s)
	for (auto i = 1; i <= n; ++i)
	{
		const auto part1 = pow(-1, i - 1);
		const auto part3 = (2 * i) - 1;
		const auto part2 = pow(x, part3);
		// cout << omp_get_thread_num() <<endl;
		s +=  part1 * (part2 / part3);		
	}

	return s;
}
