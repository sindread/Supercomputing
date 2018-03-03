#pragma once

const double mach1_expected_value_after_3_iterations = 3.141621029325034625046832517116408069706244618213430554860; //3.141621029325034425046832517116408069706244618213430554860

double mach1_calculate_pi(const int &rank, const int &numberOfProcessors, const int &numberOfIntervals);
double arctan(const int& numberOfProcessors, const int& n, const double& x, const int& rank, double* answers);
