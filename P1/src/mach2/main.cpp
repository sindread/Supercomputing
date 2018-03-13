#include "mach2.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <limits>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

#define M_PI acos(-1.0)

using namespace std;

int main(int argc, char* argv[])
{
	auto maxK = 7;
	auto plot = false;

	if (argc > 1){
		string arg = argv[1];
		if (arg =="-v"){
			plot = true;
		}
	} 

	if (plot) {
		for (int k = 1; k <= maxK ; k++){
			auto n = pow (2, k);
			auto start = Clock::now(); 

			auto pi = mach2_calculate_pi(n);

			auto end = Clock::now(); 
			
			cout << "Pi is with mach2, with " << n << " iterations: Pi_" << n <<" = "<< pi <<  endl;
			cout << "Error(PI-pi_" << n << "): E  = " << abs(M_PI-pi) <<  endl;
			cout << "Runtime: Time = " << chrono::duration_cast<chrono::nanoseconds>(end - start).count() << "ns" << endl;
		}
	}
	else {
		auto n = 1000;
		auto start = Clock::now(); 

		auto pi = mach2_calculate_pi(n);

		auto end = Clock::now(); 
		
		cout << "Pi is with mach2, with " << n << " iterations: Pi_" << n <<" = "<< pi <<  endl;
		cout << "Error(PI-pi_" << n << "): E  = " << abs(M_PI-pi) <<  endl;
		cout << "Runtime: Time = " << chrono::duration_cast<chrono::nanoseconds>(end - start).count() << "ns" << endl;
	}
	
		

	return 0;
}