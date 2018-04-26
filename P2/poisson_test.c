#include "poisson_test.h"

bool assert_equal(int a, int b){
    if (a != b){
        return false;
    }
}

void unit_transpose_parallel(){
    real** arr = mk_2D_array(3,3,false);
    printMatrix(arr, 3);
}

void lengthOfWork_FourWorkersTenWorkTasks_AllResiveTheSame(int m, int numProcs, int rank){
    int *sendcounts, *sdispls, *recvcounts, *rdispls;
    int expectedSendcounts[] = {3*10, 3*10, 2*10, 2*10}; 
    int expectedSdispls[] = {0*10, 3*10, 6*10, 8*10};
    int expectedRecvcounts[] = {3, 3, 2, 2};
    int expectedRdispls[]  = {0, 3, 6, 8};

    length_of_work(m, numProcs, rank);

    printf ("Unit test for processor: %d", rank);
    
    for (int i = 0; i < numProcs ; i++){
        if (!assert_equal(expectedSendcounts[rank], sendcounts[i])){
            printf("Processor %d : Not expected sendcounts in length of work, expeced: %d was %d", rank , expectedSendcounts[rank], sendcounts[i]);
        }   
        if (!assert_equal(expectedSdispls[rank], sdispls[i])){
            printf("Processor %d : Not expected sendcounts in length of work, expeced: %d was %d", rank , expectedSdispls[rank], sdispls[i]);
        }  
        if (!assert_equal(expectedRecvcounts[i], recvcounts[i])){
            printf("Processor %d : Not expected recvcounts in length of work, expeced: %d was %d", rank , expectedRecvcounts[i], recvcounts[i]);
        } 
        if (!assert_equal(expectedRdispls[i], rdispls[i])){
            printf("Processor %d : Not expected rdispls in length of work, expeced: %d was %d", rank , expectedRdispls[i], rdispls[i]);
        } 
    }
}

void run_poisson_unit_tests(int numProcs, int rank, int m){
    unit_transpose_parallel();

    lengthOfWork_FourWorkersTenWorkTasks_AllResiveTheSame(m, numProcs, rank);
}
