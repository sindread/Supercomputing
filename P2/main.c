
#include <mpi.h>
#include "poisson.h"
#include "poisson_test.h"

int main(int argc, char **argv) {
    int numProcs, rank, numThreads, n;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Min rank %d \n" , rank);
    
    n = atoi(argv[1]);
    numThreads = atoi(argv[2]);

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
        
        printf("Running with %d processes %d threads\n", numProcs, numThreads);
    }

    if (argc == 3){
        if(rank == 0){
            printf("Kjører program \n");
        }

        run_poisson(numProcs, rank, numThreads, n);
    }
    else{
        if(rank == 0){
            printf("Kommer i test \n");
        }
        //run_poisson_unit_tests(numProcs, rank, n-1);
    }

    MPI_Finalize();
    return 0;
}