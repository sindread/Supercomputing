#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define true 1
#define false 0

typedef double real;
typedef int bool;

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
void transpose_parallel(real **bt, real **b, size_t m, int *sendcounts, int *sdispls, int *recvcounts, int *rdispls, MPI_Datatype mpi_matrix);
void printMatrix(real** matrix, int size);
void create_mpi_datatype(size_t m, MPI_Datatype mpi_vector, MPI_Datatype mpi_matrix);
void free_mpi_datatype(MPI_Datatype mpi_vector, MPI_Datatype mpi_matrix);
void length_of_work(int m, int numProcs, int rank, int* sendcounts, int* sdispls, int* recvcounts, int* rdispls);
int run_poisson(int numProcs, int rank, int numThreads, int n);