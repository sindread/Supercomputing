/**
 * C program to solve the two-dimensional Poisson equation on
 * a unit square using one-dimensional eigenvalue decompositions
 * and fast sine transforms.
 *
 * Einar M. Rønquist
 * NTNU, October 2000
 * Revised, October 2001
 * Revised by Eivind Fonn, February 2015
 */

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

// Function prototypes
real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
void transpose(real **bt, real **b, size_t m);
real rhs(real x, real y);

// Functions implemented in FORTRAN in fst.f and called from C.
// The trailing underscore comes from a convention for symbol names, called name
// mangling: if can differ with compilers.
void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);


//Our code
void transpose_parallel(real **bt, real **b, size_t m, int* sendcount, int* senddis, int* recvcount, int* recvdis);
void divide_work(size_t m, int numProcs, int rank, int* sendr, int* senddis, int* recvr, int* recvdis);

int main(int argc, char **argv)
{
    int numProcs, rank, numThreads, startTime;
    int* sendcount, *senddis, *recvcount, *recvdis;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    numThreads = atoi(argv[2]);

    printf("numprocs %d \n", numProcs);
    if(rank == 0){
        if (argc < 3) {
            printf("Usage:\n");
            printf("  poisson n\n\n");
            printf("Arguments:\n");
            printf("  n: the problem size (must be a power of 2)\n");
            printf("  t: number of threads\n");

            MPI_Finalize();
            return -1;
        }
        
        printf("Running with %d processes and %d threads\n", numProcs, numThreads);
        startTime = MPI_Wtime();
    }
    
    
    omp_set_dynamic(0);
    omp_set_num_threads(numThreads);
    

    /*
     *  The equation is solved on a 2D structured grid and homogeneous Dirichlet
     *  conditions are applied on the boundary:
     *  - the number of grid points in each direction is n+1,
     *  - the number of degrees of freedom in each direction is m = n-1,
     *  - the mesh size is constant h = 1/n.
     */
    int n = atoi(argv[1]);
    int m = n - 1;
    real h = 1.0 / n;

    if((n & (n-1))!= 0 && n){
        if(rank == 0){
            printf("Problem size needs to be power of 2\n");
        }

        MPI_Finalize();
        return -1;
    }
    sendcount = (int*) malloc(numProcs*sizeof(int));
    recvcount = (int*) calloc(numProcs,sizeof(int));
    senddis = (int*) malloc(numProcs*sizeof(int));
    recvdis = (int*) calloc(numProcs,sizeof(int));

    divide_work(m, numProcs, rank, sendcount, senddis, recvcount, recvdis);

    /*
     * Grid points are generated with constant mesh size on both x- and y-axis.
     */
    real *grid = mk_1D_array(n+1, false);
    #pragma omp parallel for 
    for (size_t i = 0; i < n+1; i++) {
        grid[i] = i * h;
    }

    /*
     * The diagonal of the eigenvalue matrix of T is set with the eigenvalues
     * defined Chapter 9. page 93 of the Lecture Notes.
     * Note that the indexing starts from zero here, thus i+1.
     */
    real *diag = mk_1D_array(m, false);
    #pragma omp parallel for 
    for (size_t i = 0; i < m; i++) {
        diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
    }

    /*
     * Allocate the matrices b and bt which will be used for storing value of
     * G, \tilde G^T, \tilde U^T, U as described in Chapter 9. page 101.
     */
    real **b = mk_2D_array(m, m, false);
    real **bt = mk_2D_array(m, m, false);

    /*
     * This vector will holds coefficients of the Discrete Sine Transform (DST)
     * but also of the Fast Fourier Transform used in the FORTRAN code.
     * The storage size is set to nn = 4 * n, look at Chapter 9. pages 98-100:
     * - Fourier coefficients are complex so storage is used for the real part
     *   and the imaginary part.
     * - Fourier coefficients are defined for j = [[ - (n-1), + (n-1) ]] while 
     *   DST coefficients are defined for j [[ 0, n-1 ]].
     * As explained in the Lecture notes coefficients for positive j are stored
     * first.
     * The array is allocated once and passed as arguments to avoid doings 
     * reallocations at each function call.
     */
    int nn = 4 * n;
    real *z = mk_1D_array(nn, false);

    /*
     * Initialize the right hand side data for a given rhs function.
     * Note that the right hand-side is set at nodes corresponding to degrees
     * of freedom, so it excludes the boundary (bug fixed by petterjf 2017).
     * 
     */
    #pragma omp parallel for  collapse(2)
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            b[i][j] = h * h * rhs(grid[i+1], grid[j+1]);
        }
    }

    /*
     * Compute \tilde G^T = S^-1 * (S * G)^T (Chapter 9. page 101 step 1)
     * Instead of using two matrix-matrix products the Discrete Sine Transform
     * (DST) is used.
     * The DST code is implemented in FORTRAN in fsf.f and can be called from C.
     * The array zz is used as storage for DST coefficients and internally for 
     * FFT coefficients in fst_ and fstinv_.
     * In functions fst_ and fst_inv_ coefficients are written back to the input 
     * array (first argument) so that the initial values are overwritten.
     */
    #pragma omp parallel for 
    for (size_t i = 0; i < m; i++) {
        fst_(b[i], &n, z, &nn);
    }
    //transpose(bt, b, m); Need to implement transpose for parallel
    transpose_parallel(bt, b, m, sendcount, senddis, recvcount, recvdis);
    #pragma omp parallel for 
    for (size_t i = 0; i < m; i++) {
        fstinv_(bt[i], &n, z, &nn);
    }

    /*
     * Solve Lambda * \tilde U = \tilde G (Chapter 9. page 101 step 2)
     */
    #pragma omp parallel for 
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            bt[i][j] = bt[i][j] / (diag[i] + diag[j]);
        }
    }

    /*
     * Compute U = S^-1 * (S * Utilde^T) (Chapter 9. page 101 step 3)
     */
    #pragma omp parallel for 
    for (size_t i = 0; i < m; i++) {
        fst_(bt[i], &n, z, &nn);
    }
    //transpose(b, bt, m);
    transpose_parallel(bt, b, m, sendcount, senddis, recvcount, recvdis);
    
    #pragma omp parallel for 
    for (size_t i = 0; i < m; i++) {
        fstinv_(b[i], &n, z, &nn);
    }

    /*
     * Compute maximal value of solution for convergence analysis in L_\infty
     * norm.
     */
    double u_max = 0.0;
    #pragma omp parallel for  collapse(2)
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            u_max = u_max > b[i][j] ? u_max : b[i][j];
        }
    }
    
    if(rank == 0){
        printf("u_max = %e\n", u_max);    
    }
    return 0;
}

/*
 * This function is used for initializing the right-hand side of the equation.
 * Other functions can be defined to swtich between problem definitions.
 */

real rhs(real x, real y) 
{
    return 2 * (y - y*y + x - x*x);
}

/*
 * Write the transpose of b a matrix of R^(m*m) in bt.
 * In parallel the function MPI_Alltoallv is used to map directly the entries
 * stored in the array to the block structure, using displacement arrays.
 */

void transpose(real **bt, real **b, size_t m)
{
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            bt[i][j] = b[j][i];
        }
    }
}

/*
 * The allocation of a vectore of size n is done with just allocating an array.
 * The only thing to notice here is the use of calloc to zero the array.
 */

real *mk_1D_array(size_t n, bool zero)
{
    if (zero) {
        return (real *)calloc(n, sizeof(real));
    }
    return (real *)malloc(n * sizeof(real));
}

/*
 * The allocation of the two-dimensional array used for storing matrices is done
 * in the following way for a matrix in R^(n1*n2):
 * 1. an array of pointers is allocated, one pointer for each row,
 * 2. a 'flat' array of size n1*n2 is allocated to ensure that the memory space
 *   is contigusous,
 * 3. pointers are set for each row to the address of first element.
 */

real **mk_2D_array(size_t n1, size_t n2, bool zero)
{
    // 1
    real **ret = (real **)malloc(n1 * sizeof(real *));

    // 2
    if (zero) {
        ret[0] = (real *)calloc(n1 * n2, sizeof(real));
    }
    else {
        ret[0] = (real *)malloc(n1 * n2 * sizeof(real));
    }
    
    // 3
    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }
    return ret;
}

/*
 *
 *  OUR CODE
 * 
 */

void transpose_parallel(real **bt, real **b, size_t m, int* sendcount, int* senddis, int* recvcount, int* recvdis){
    MPI_Alltoallv(b[0], sendcount, senddis, MPI_DOUBLE, bt[0], recvcount, recvdis, MPI_DOUBLE, MPI_COMM_WORLD);
}



void divide_work(size_t m, int numProcs, int rank, int* sendcount, int* senddis, int* recvcount, int* recvdis){
    int rows_pr_processor = m/numProcs;
    int remainder = m%numProcs;

    for(int i=0; i<numProcs; i++){
        if(m <= numProcs && i < m){
            recvcount[i] = 1;
            recvdis[i]=i;
        }else if(remainder != 0){
            if(remainder > i){
                recvcount[i] = rows_pr_processor + 1;
                recvdis[i] = i*rows_pr_processor + i;
            }else{
                recvcount[i] = rows_pr_processor;
                recvdis[i] = i*rows_pr_processor + remainder;
            }
        }
    }

    for(int i=0; i<numProcs; i++){
        sendcount[i]=recvcount[rank]*m;
        senddis[i]=recvcount[rank]*m;
    }
}