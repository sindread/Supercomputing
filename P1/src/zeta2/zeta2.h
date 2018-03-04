#pragma once

const double zeta2_expected_value_after_3_iterations = 2.85773803324704145; //114563498087156873290626938727432781816504;

double zeta2_calculate_pi(const int &rank, const int &numberOfProcessors, const int &numberOfIntervals);
void zeta2_calculate_vi(const int &numberOfProcesses, const int &rank, const int &n, double &sumPiParts);

