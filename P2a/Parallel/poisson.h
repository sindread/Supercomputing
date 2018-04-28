#ifndef POISSON_H
#define POISSON_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define true 1
#define false 0

typedef double real;
typedef int bool;

int *sendcounts, *sdispls, *recvcounts, *rdispls;

// Functions implemented in FORTRAN in fst.f and called from C.
// The trailing underscore comes from a convention for symbol names, called name
// mangling: if can differ with compilers.
void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);

real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
void transpose(real **bt, real **b, size_t m);
real rhs(real x, real y);

void printMatrix(real** matrix, int length);
int* mk_1D_array_int(size_t n, bool zero);
void transpose_parallel(real **bt, real **b, size_t m);
void printMatrix(real** matrix, int size);
void create_mpi_datatype(size_t m);
void free_mpi_datatype();
void transpose_parallel_setup(int m, int numProcs, int rank);
void run_poisson(int numProcs, int rank, int numThreads, int n);

#endif