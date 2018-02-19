#include "stdafx.h"
#include "zeta0.h"
#include <math.h>

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

const bool zeta0::zUnitTest()
{
	double answerShouldBe = 3.082207001484488225125096190727122112617812011722287272;

	return answerShouldBe == calculatePi(3);
}

