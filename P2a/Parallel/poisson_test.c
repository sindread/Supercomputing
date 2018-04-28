#include "poisson_test.h"

int unit_transpose_parallel(int numProcs, int rank){
    int m = 3;

    real** matrix = mk_2D_array(m, m, true);
    real** matrix_transposed = mk_2D_array(m, m, true);
    real** expected_matrix_tranpose = mk_2D_array(m, m, true);

    int number = 1;

    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            matrix[i][j] = number;
            expected_matrix_tranpose [j][i] = number;
            number++;
        }
    }

    transpose_parallel_setup(m, numProcs, rank);
    transpose_parallel(matrix, matrix_transposed, m);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if (matrix_transposed[i][j] - expected_matrix_tranpose[i][j] != 0){
                return -1;
            }
        }   
    }     

    return 1;
}

int lengthOfWork_FourWorkersTenWorkTasks_AllResiveTheSame(int m, int numProcs, int rank){
    int expectedSendcounts[] = {50, 50}; 
    int expectedSdispls[] = {0, 50};
    int expectedRecvcounts[] = {5, 5};
    int expectedRdispls[]  = {0, 5};

    transpose_parallel_setup(m, numProcs, rank);

    //printf ("Unit test for processor: %d \n", rank);
    
    for (int i = 0; i < numProcs ; i++){
        if (expectedSendcounts[rank] - sendcounts[i] != 0){
            //printf("Processor %d : Not expected sendcounts in length of work, expeced: %d was %d \n", rank , expectedSendcounts[rank], sendcounts[i]);
            return -1;
        }   
        if (expectedSdispls[rank] - sdispls[i] != 0){
            //printf("Processor %d : Not expected sendcounts in length of work, expeced: %d was %d \n", rank , expectedSdispls[rank], sdispls[i]);
            return -1;
        }  
        if (expectedRecvcounts[i] - recvcounts[i] != 0){
            //printf("Processor %d : Not expected recvcounts in length of work, expeced: %d was %d \n", rank , expectedRecvcounts[i], recvcounts[i]);
            return -1;
        } 
        if (expectedRdispls[i] - rdispls[i] != 0){
            //printf("Processor %d : Not expected rdispls in length of work, expeced: %d was %d \n", rank , expectedRdispls[i], rdispls[i]);
            return -1;
        } 
    }
    return 1;
}

void run_poisson_unit_tests(int numProcs, int rank, int m){
    int failed = 0;
    int passed = 0;

    if (unit_transpose_parallel(numProcs, rank) > 0){
        printf("unit_transpose_parallel for process %d: Pass\n", rank);
        passed++;
    }else {
        printf("unit_transpose_parallel for process %d: Failed\n", rank);
        failed++;
    }

    if (lengthOfWork_FourWorkersTenWorkTasks_AllResiveTheSame(m, numProcs, rank) > 0){
        printf("lengthOfWork_FourWorkersTenWorkTasks_AllResiveTheSame for process %d: Pass\n", rank);
        passed++;
    } else {
        printf("lengthOfWork_FourWorkersTenWorkTasks_AllResiveTheSame for process %d: Failed\n", rank);
        failed++;
    }

    printf("Tests passes %d of %d for process %d", passed, passed + failed, rank);
}
