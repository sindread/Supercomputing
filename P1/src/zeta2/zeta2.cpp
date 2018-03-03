#include "zeta2.h"
#include <cmath>
#include <omp.h>
#include <chrono>
#include <iostream>

using namespace std;

void zeta2_calculate_vi(const int &n, double &answer)
{

	#pragma omp parallel for reduction (+:answer)
	for (auto i = 1; i <= n; ++i)
	{	
		cout << omp_get_thread_num() <<endl;
		answer += 1.0 / pow(i,2);
	}
	
	
}


double zeta2_calculate_pi(const int &numberOfIntervals){
	
	double pi = 0.0;
	auto startTime = chrono::high_resolution_clock::now();

	zeta2_calculate_vi(numberOfIntervals, pi);

	pi = sqrt(6*pi);

	auto stopTime = chrono::high_resolution_clock::now();


	return pi;
	
}



