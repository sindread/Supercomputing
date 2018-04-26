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

#include "poisson.h"

void run_poisson(int numProcs, int rank, int numThreads, int n){
    /*
     *  The equation is solved on a 2D structured grid and homogeneous Dirichlet
     *  conditions are applied on the boundary:
     *  - the number of grid points in each direction is n+1,
     *  - the number of degrees of freedom in each direction is m = n-1,
     *  - the mesh size is constant h = 1/n.
     */

    omp_set_dynamic(0);
    omp_set_num_threads(numThreads);
    
    int m = n - 1;
    real h = 1.0 / n;

    if((n & (n-1))!= 0 && n){
        if(rank == 0){
            printf("Problem size needs to be power of 2\n");
        }

        MPI_Finalize();
    }

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

    length_of_work(m, numProcs, rank);

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

    create_mpi_datatype(m);

    //transpose(bt, b, m); Need to implement transpose for parallel
    transpose_parallel(b, bt, m);

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

    if(rank == 0){
        printf("test 1 \n");    
    }

    /*
     * Compute U = S^-1 * (S * Utilde^T) (Chapter 9. page 101 step 3)
     */
    #pragma omp parallel for 
    for (size_t i = 0; i < m; i++) {
        fst_(bt[i], &n, z, &nn);
    }
    //transpose(b, bt, m);
    transpose_parallel(bt, b, m);
    
    #pragma omp parallel for 
    for (size_t i = 0; i < m; i++) {
        fstinv_(b[i], &n, z, &nn);
    }
    if(rank == 0){
        printf("test 2 \n");    
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
        printMatrix(b,m);
    }

    free_mpi_datatype();
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

int *mk_1D_array_int(size_t n, bool zero)
{
    if (zero) {
        return (int *)calloc(n, sizeof(int));
    }
    return (int *)malloc(n * sizeof(int));
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
MPI_Datatype mpi_vector;
MPI_Datatype mpi_matrix;

void create_mpi_datatype(size_t m){
    //Creating datatype for vectors for mpi
    MPI_Type_vector(m, 1, m, MPI_DOUBLE , &mpi_vector);
    MPI_Type_commit(&mpi_vector);

    //Using vectors as columns in matrix
    MPI_Type_create_resized(mpi_vector, 0, sizeof(double),&mpi_matrix);
    MPI_Type_commit(&mpi_matrix);
}

void free_mpi_datatype(){
    MPI_Type_free(&mpi_vector);
    MPI_Type_free(&mpi_matrix);
}

void transpose_parallel(real **b, real **bt, size_t m){
    MPI_Alltoallv(b[0], sendcounts, sdispls, MPI_DOUBLE, bt[0], recvcounts, rdispls, mpi_matrix, MPI_COMM_WORLD);
}

void length_of_work(int m, int numProcs, int rank){
    
    sendcounts = mk_1D_array_int(numProcs, false);
    sdispls = mk_1D_array_int(numProcs, false);
    recvcounts = mk_1D_array_int(numProcs, false);
    rdispls = mk_1D_array_int(numProcs, true);

    int workTasks = m/numProcs;
    int workLeft = m%numProcs;

    for (int i = 0; i < numProcs; i++){
        recvcounts[i] = workTasks;
        
        if (workLeft > 0){
            recvcounts[i]++;
            workLeft--;
        }  

        if (i > 0){ 
            rdispls[i] = rdispls[i-1] + recvcounts[i-1];
        }
    }

    for (int j = 0; j<numProcs; j++){
        sendcounts[j] = recvcounts[rank]*m;
        sdispls[j] = rdispls[rank]*m; 
    }
}

void printMatrix(real** matrix, int length){
    for(size_t i = 0; i < length; i++){
        for(size_t j = 0; j < length; j++){
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}