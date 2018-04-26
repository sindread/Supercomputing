#pragma once
#include <stdlib.h>
#include <stdio.h>
#include "poisson.h"

#define true 1
#define false 0

typedef double real;
typedef int bool;

void run_poisson_unit_tests(int numProcs, int rank, int m);