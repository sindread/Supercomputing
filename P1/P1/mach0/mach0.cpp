#include "stdafx.h"
#include "mach0.h"
#include <iostream>

#define M_PI acos(-1.0)

const double mach0:: calculatePi(const int &n)
{
	const double x1 = 1.0 / 5;
	const double x2 = 1.0 / 239;

	double s1 = 0, s2 = 0;

	for (int i = 1; i <= n; ++i)
	{
		s1 += calcV(i, x1);
		s2 += calcV(i, x2);
	}

	return 4 * (4 * s1 - s2);
}


const bool mach0::mUnitTest()
{
	const double answerShouldBe = 3.141621029325034625046832517116408069706244618213430554860; //3.141621029325034425046832517116408069706244618213430554860

	return answerShouldBe == calculatePi(3);
}


const double mach0::calcV(const double &i, const double &x)
{
	const int part1 = pow(-1, i-1);
	const int part3 = 2 * i - 1;
	const double part2 = pow(x, part3);
	return part1 * (part2 / part3);
}


double* mach0::mVerificationTest()
{
	double answers[24];
	for (int i = 1; i <= 24; ++i)
	{
		const int n = pow(2, i);
		answers[i - 1] = abs(M_PI - calculatePi(n));
	}

	return answers;
}



