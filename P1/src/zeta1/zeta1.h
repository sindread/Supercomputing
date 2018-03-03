#pragma once

const double zeta1_expected_value_after_3_iterations = 2.85773803324704145; //114563498087156873290626938727432781816504;

double zeta1_calculate_pi(const int &rank, const int &numberOfProcessors, const int &numberOfIntervals);
void zeta1_calculate_vi(const int &numberOfProcesses, const int &rank, const int &n, double &sumPiParts);
void sumVector(const double* vector, const int& length, double& sum);
