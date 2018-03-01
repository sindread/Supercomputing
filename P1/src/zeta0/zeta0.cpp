#include "zeta0.h"
#include <cmath>

double zeta0_calculate_pi(const int &n)
{
	auto s = 0.0;

	for (auto i = 1; i <= n; ++i)
	{
		s += 1.0 / pow(i,2);
	}

	//pi = (6*s)^(1/2);
	return sqrt(s * 6);
}