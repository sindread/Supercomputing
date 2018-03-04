#pragma once

const double zeta1_expected_value_after_3_iterations = 2.85773803324704145; //114563498087156873290626938727432781816504;
const int TAG_LENGTH = 42;
const int TAG_VPARTS = 1337;

void master_init(int argc, char* argv[], int &n, bool &isUnitTest);
void master_task(const int &n, const bool &isUnitTest, const int &numberOfProcesses);
void slave_task(int &rank, int &numberOfProcesses);
void vi_parts(const int &n, double* vi);
void sumVector(const double* vector, const int& length, double& sum);