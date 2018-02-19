#include "stdafx.h"
#include "zeta0.h"
#include <iostream>

const double zeta0::calculatePi(const int &n)
{
	double s = 0.0;

	for (int i = 1; i <= n; ++i)
	{
		const double v_i = 1. / (i*i);
		s += v_i;
	}

	return sqrt(s * 6);
}
