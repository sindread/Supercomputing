// P1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "zeta0.h"
#include "mach0.h"
#include <iostream>
#include <limits>

using namespace std;


int main()
{
	zeta0 zeta;
	mach0 mach;
	
	cout.precision(std::numeric_limits<double>::digits10 + 2);
	cout << fixed << zeta.calculatePi(10000) << endl;
	cout << fixed << mach.calculatePi(10000) << endl;

	cout << zeta.zUnitTest() << endl;

    return 0;
}

