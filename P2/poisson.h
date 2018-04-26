#include <stdio.h>

typedef double real;
typedef int bool;

real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
void transpose(real **bt, real **b, size_t m);
real rhs(real x, real y);

void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);

void transpose_parallel(real **bt, real **b, size_t m);
void divide_work(size_t m, int numProcs, int rank);
void printMatrix(real** matrix, int size);
void create_mpi_datatype(size_t m);
void free_mpi_datatype();