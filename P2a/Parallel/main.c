
#include <mpi.h>
#include <string.h>
#include "poisson.h"
#include "poisson_test.h"

void validate(int rank, int numProcs);

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
        
        //printf("Running with %d processes %d threads:\n", numProcs, numThreads);
    }

    if (argc == 3){
        if(rank == 0){
            //(printf("Running program \n");
            startTime = MPI_Wtime();
        }
        
        MPI_Barrier(MPI_COMM_WORLD);

        run_poisson(numProcs, rank, numThreads, n);

        if(rank == 0){
            double endTime = MPI_Wtime() - startTime;
            printf("Elapsed time: %f ms\n", endTime*1000);
        }
    }
    else{
        if (strcmp (argv[3],"u") != 1){
            if (numProcs == 2) {
                if(rank == 0 ){
                    printf("Running unit-tests: \n");
                    startTime = MPI_Wtime();
                }
                MPI_Barrier(MPI_COMM_WORLD);

                run_poisson_unit_tests(numProcs, rank, 10);
            } 
            else {
                if(rank == 0 ){
                    printf("Running unit-tests must have 2 processes, was %d \n", numProcs);
                }
                
            }
        }else {
            validate(rank, numProcs);
        }
    }
    
    MPI_Finalize();
    return 0;
}

void validate(int rank, int numProcs){
    if (rank == 0) {
        printf("Running validation-test for %d processes: \n", numProcs);
    }

    for (int t = 0; t < 6+1; t++){
        for (int k = 1; k < 14 +1; k++){
            double startTime = MPI_Wtime();
            int nt = pow(2, t);
            int n = pow (2, k);
            //double u_max = 
            run_poisson(numProcs, rank, nt, n);

            if(rank == 0){
                double endTime = MPI_Wtime() - startTime;
                printf("%d threads with problem size %d, elapsed time: %f ms\n", nt, n, endTime*1000);
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}