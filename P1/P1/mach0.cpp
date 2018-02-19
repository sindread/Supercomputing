#include "stdafx.h"
#include "mach0.h"
#include <math.h>


const double mach0:: calculatePi(const int &n)
{
	const double x1 = 1.0 / 5;
	const double x2 = 1.0 / 239;

	double s1 = 0, s2 = 0;

	for (int i = 1; i < n; ++i)
	{
		s1 += calcV(i, x1);
		s2 += calcV(i, x2);
	}

	return 4 * (4 * s1 - s2);
}


const double mach0::calcV(const double &i, const double &x)
{
	const int part1 = pow(-1, i-1);
	const int part3 = 2 * i - 1;
	const double part2 = pow(x, part3);
	return part1 * (part2 / part3);
}
