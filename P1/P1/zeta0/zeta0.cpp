#include "stdafx.h"
#include "zeta0.h"
#include <iostream>

#define M_PI acos(-1.0)

const double zeta0::calculatePi(const int &n)
{
	double s = 0.0;

	for (int i = 1; i <= n; ++i)
	{
		const double v_i = 1.0 / (i*i);	
		s += v_i;
	}

	return sqrt(s * 6);
}

const bool zeta0::zUnitTest()
{
	const double answerShouldBe = 2.85773803324704145; //114563498087156873290626938727432781816504;

	return answerShouldBe == calculatePi(3);
}


double* zeta0::zVerificationTest()
{
	double answers[24];
	for (int i = 1; i <= 24; ++i)
	{
		const int n = pow(2, i);
		std::cout << abs(M_PI - calculatePi(n)) << std::endl;
		answers[i - 1] = abs(M_PI -  calculatePi(n));
	}

	return answers;
}


