#pragma once

const double mach2_expected_value_after_3_iterations = 3.141621029325034625046832517116408069706244618213430554860; //3.141621029325034425046832517116408069706244618213430554860

double mach2_calculate_pi(const int &rank, const int &numberOfProcessors, const int &numberOfIntervals);
void mach2_calculate_vi(const int& numberOfProcessors, const int& n, const int& rank, double& answer);
double arctan(const int& numberOfProcessors, const int& n, const double& x, const int& rank);
