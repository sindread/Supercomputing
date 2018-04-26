
#include <mpi.h>
#include "poisson.h"
#include "poisson_test.h"

int main(int argc, char **argv) {
    int numProcs, rank, numThreads, n;
    double startTime = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
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
        
        printf("Running with %d processes %d threads:\n", numProcs, numThreads);
    }

    if (argc == 3){
        if(rank == 0){
            printf("Running program \n");
            startTime = MPI_Wtime();
        }
        
        MPI_Barrier(MPI_COMM_WORLD);

        run_poisson(numProcs, rank, numThreads, n);
    }
    else{
        if (numProcs == 2) {
            if(rank == 0 ){
                printf("Running unit-tests: \n");
            }
            MPI_Barrier(MPI_COMM_WORLD);

            run_poisson_unit_tests(numProcs, rank, 10);
        } 
        else {
            if(rank == 0 ){
                printf("Running unit-tests must have 2 processes, was %d \n", numProcs);
            }
            
            MPI_Finalize();
            return -1;
        }
    }

    if(rank == 0){
        double endTime = MPI_Wtime() - startTime;
        printf("Elapsed time: %f ms\n", endTime*1000);
    }
    MPI_Finalize();
    return 0;
}