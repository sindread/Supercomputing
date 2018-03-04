#include "mach0.h"
#include <cmath>

double mach0_calculate_pi(const int& n) 
{
	const auto x1 = double(1)/5;
	const auto x2 = double(1)/239;

	// pi = 4 * (4 * arctan(1/5) - arctan (1/239))
	return 4 * (4 * arctan(n, x1) - arctan(n, x2));
}

double arctan(const int& n, const double& x)
{
	auto s = 0.0;

	for (auto i = 1; i <= n; ++i)
	{
		const auto part1 = pow(-1, i - 1);
		const auto part3 = (2 * i) - 1;
		const auto part2 = pow(x, part3);
		s += part1 * (part2 / part3);
	}

	return s;
}


