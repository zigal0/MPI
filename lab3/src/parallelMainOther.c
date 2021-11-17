//
// Created by zigal0 on 16.11.2021.
//
// Other realization parallel main task (master calculate as well + all in main)
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

// size of array to calculate
#define ISIZE 1000
#define JSIZE 1000

// tags for send/receive functions
#define TAG_START 1
#define TAG_SIZE 2
#define TAG_RESULT 3

void print_to_file(FILE *ff, double **a);

int main(int argc, char **argv) {
    // variables
    int rank = 0, size = 0;         // rank of process and quantity of processes
    int i, j;                       // for iteration
    MPI_Status status;              // status for checking delivery
    MPI_Request request;            // status of delivery
    // MPI initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // master
    if (rank == 0) {
        // variables
        double time_start, time_finish;                     // for measuring time
        int curNum = ISIZE;                                 // current number of columns
        double **a = malloc(ISIZE * sizeof(double *));  // array to calculate
        for (i = 0; i < ISIZE; i++) {
            a[i] = malloc(JSIZE * sizeof(double));
        }
        int *localSize = (int *) malloc(size * sizeof(int));  // for uniform distribution between processes
        int *localPos = (int *) malloc(size * sizeof(int));   // start position for iteration in every process
        localPos[0] = 0;
        time_start = MPI_Wtime();
        // distribution
        for (i = 0; i < size; ++i) {
            if (curNum % (size - i) != 0) {
                localSize[i] = curNum / (size - i) + 1;
            } else {
                localSize[i] = curNum / (size - i);
            }
            curNum -= localSize[i];
            if (i != 0) {
                MPI_Isend(&localSize[i], 1, MPI_INT, i, TAG_SIZE, MPI_COMM_WORLD, &request);
                MPI_Isend(&localPos[i], 1, MPI_INT, i, TAG_START, MPI_COMM_WORLD, &request);
            }
            if (i != size - 1) {
                localPos[i + 1] = localPos[i] + localSize[i];
            }
        }
        // own calculation
        for (i = localPos[0]; i < localSize[0]; ++i) {
            for (j = 0; j < JSIZE; ++j) {
                a[i][j] = sin(0.00001 * (10 * i + j));
            }
        }
        // collect all parts of task
        for (i = 1; i < size; ++i) {
            for (j = localPos[i]; j < localPos[i] + localSize[i]; ++j) {
                MPI_Recv(a[j], JSIZE, MPI_DOUBLE, i, TAG_RESULT, MPI_COMM_WORLD, &status);
            }
        }
        time_finish = MPI_Wtime(); // time of finish
        printf("Execution time - %0.6f\n", time_finish - time_start);


        // print data to txt
        FILE *ff = fopen("output/resultPMO.txt", "w");
        print_to_file(ff, a);
        fclose(ff);

        // free all memory
        free(localSize);
        free(localPos);
        for (i = 0; i < ISIZE; i++) {
            free(a[i]);
        }
        free(a);

        // slaves
    } else {
        //variables
        int localSize = 0;      // size of batch
        int localPos = 0;       // start index in general array
        MPI_Recv(&localSize, 1, MPI_INT, 0, TAG_SIZE, MPI_COMM_WORLD, &status);
        MPI_Recv(&localPos, 1, MPI_INT, 0, TAG_START, MPI_COMM_WORLD, &status);
        // memory for result
        double **res = (double **) malloc(localSize * sizeof(double *));
        for (i = 0; i < localSize; ++i) {
            res[i] = (double *) malloc(JSIZE * sizeof(double));
        }
        // calculation
        for (i = 0; i < localSize; ++i) {
            for (j = 0; j < JSIZE; ++j) {
                res[i][j] = sin(0.00001 * (10 * (i + localPos) + j));
            }
        }
        // send part of task to master
        for (i = 0; i < localSize; ++i) {
            MPI_Send(res[i], JSIZE, MPI_DOUBLE, 0, TAG_RESULT, MPI_COMM_WORLD);
        }
        // free all structures except res
        for (i = 0; i < localSize; ++i) {
            free(res[i]);
        }
        free(res);
    }

    // waiting for all processes
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

void print_to_file(FILE *ff, double **a) {
    for (int i = 0; i < ISIZE; i++) {
        for (int j = 0; j < JSIZE; j++) {
            fprintf(ff, "%f ", a[i][j]);
        }
        fprintf(ff, "\n");
    }
}
