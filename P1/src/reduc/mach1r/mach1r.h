#pragma once

const double mach1_expected_value_after_3_iterations = 3.141621029325034625046832517116408069706244618213430554860; //3.141621029325034425046832517116408069706244618213430554860

const int TAG_LENGTH = 42;
const int TAG_VPARTS_5 = 1337;
const int TAG_VPARTS_239 = 13;
const int TAG_PARTSUM1 = 666;
const int TAG_PARTSUM2 = 667;

void master_init(int argc, char* argv[], int &n);
void master_task(const int &n, const int &numberOfProcesses);
void slave_task(int &rank, int &numberOfProcesses);
void length_of_work(int* lengthForRank, const int &n, const int &numberOfProcesses);
void vi_parts(const int& n, double* vi5, double* vi239);
double arctan_part(const int& i, const double& x);
void sumVector(const double* vector, const int& length, double& sum);
void globalReduce(const int &numberOfProcesses, double &sum, const int *lengthForProcesses);