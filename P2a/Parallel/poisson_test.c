#include "poisson_test.h"

void unit_transpose_parallel(){
    real** arr = mk_2D_array(3,3,false);
    printMatrix(arr, 3);
}

int lengthOfWork_FourWorkersTenWorkTasks_AllResiveTheSame(int m, int numProcs, int rank){
    int expectedSendcounts[] = {50, 50}; 
    int expectedSdispls[] = {0, 50};
    int expectedRecvcounts[] = {5, 5};
    int expectedRdispls[]  = {0, 5};

    length_of_work(m, numProcs, rank);

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
    //unit_transpose_parallel();

    if (lengthOfWork_FourWorkersTenWorkTasks_AllResiveTheSame(m, numProcs, rank) > 0){
        printf("lengthOfWork_FourWorkersTenWorkTasks_AllResiveTheSame for process %d: Pass\n", rank);
    } else {
        printf("lengthOfWork_FourWorkersTenWorkTasks_AllResiveTheSame for process %d: Failed\n", rank);
    }
}
