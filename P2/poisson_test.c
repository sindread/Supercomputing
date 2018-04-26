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

    length_of_work(m, numProcs, rank, sendcounts, sdispls, recvcounts, rdispls);

    printf ("Unit test for processor:", rank);
    
    for (int i = 0; i < numProcs ; i++){
        if (!assert_equal(expectedSendcounts[rank], sendcounts[i])){
            printf("Processor", rank , ": Not expected sendcounts in length of work, expeced: ", expectedSendcounts[rank], " was ", sendcounts[i]);
        }   
        if (!assert_equal(expectedSdispls[rank], sdispls[i])){
            printf("Processor", rank , ": Not expected sdispls in length of work, expeced: ", expectedSdispls[rank], " was ", sdispls[i]);
        }  
        if (!assert_equal(expectedRecvcounts[i], recvcounts[i])){
            printf("Processor", rank , ": Not expected recvcounts in length of work, expeced: ", expectedRecvcounts[i], " was ", recvcounts[i]);
        } 
        if (!assert_equal(expectedRdispls[i], rdispls[i])){
            printf("Processor", rank , ": Not expected rdispls in length of work, expeced: ", expectedRdispls[i], " was ", rdispls[i]);
        } 
    }
}

void run_poisson_unit_tests(int numProcs, int rank, int m){
    unit_transpose_parallel();

    lengthOfWork_FourWorkersTenWorkTasks_AllResiveTheSame(m, numProcs, rank);
}
