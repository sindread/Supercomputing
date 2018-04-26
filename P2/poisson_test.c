#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include "poisson.h"

#define true 1
#define false 0

typedef double real;
typedef int bool;

bool assert_equal(int a, int b){
    if (a != b){
        return false;
    }
}


void unit_transpose_parallel(){
    real** arr = mk_2D_array(3,3,false);
    printMatrix(arr, 3);
}


void lengthOfWork_FourWorkersTenWorkTasks_AllResiveTheSame(int m, int numProcs, int rank, int *sendcounts, int *sdispls, int *recvcounts, int *rdispls){
    
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
