//
// Created by zigal0 on 16.11.2021.
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

    int rank = 0, size = 0;         // rank of process and quantity of processes
    int i, j;
    MPI_Status status; // status for checking delivery
    MPI_Request request; // status of delivery// iterations

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    // master process
    if (rank == 0) {
        double time_start, time_finish; // measure time
        double **a = malloc(ISIZE * sizeof(double *));
        for (i = 0; i < ISIZE; i++) {
            a[i] = malloc(JSIZE * sizeof(double));
        }
        int *localSize;  // for uniform distribution between processes
        localSize = (int *) malloc(size * sizeof(int));
        int *startSpace;  // start x for iteration in every process
        startSpace = (int *) malloc(size * sizeof(int));
        int curNum = ISIZE;
        startSpace[0] = 0;
        // quantity of iterations in every process and start
        time_start = MPI_Wtime(); // time of start
        for (int p = 0; p < size; ++p) {
            if (curNum % (size - p) != 0) {
                localSize[p] = curNum / (size - p) + 1;
            } else {
                localSize[p] = curNum / (size - p);
            }
            curNum -= localSize[p];
            if (p != 0) {
                MPI_Isend(&localSize[p], 1, MPI_INT, p, TAG_SIZE, MPI_COMM_WORLD, &request);
                MPI_Isend(&startSpace[p], 1, MPI_INT, p, TAG_START, MPI_COMM_WORLD, &request);
            }
            if (p != size - 1) {
                startSpace[p + 1] = startSpace[p] + localSize[p];
            }
        }
        for (i = startSpace[0]; i < localSize[0]; ++i) {
            for (j = 0; j < JSIZE; ++j) {
                a[i][j] = sin(0.00001 * (10 * i + j));
            }
        }
        // collect all parts of task
        for (i = 1; i < size; ++i) {
            for (j = startSpace[i]; j < startSpace[i] + localSize[i]; ++j) {
                MPI_Recv(a[j], JSIZE, MPI_DOUBLE, i, TAG_RESULT, MPI_COMM_WORLD, &status);
            }
        }
        time_finish = MPI_Wtime(); // time of finish
        printf("Execution time - %0.6f\n", time_finish - time_start);

//        FILE *ff = fopen("output/resultPMO10.txt", "w");
//        print_to_file(ff, a);
//        fclose(ff);

        for (i = 0; i < ISIZE; i++) {
            free(a[i]);
        }
        free(a);
    } else {
        int localSize = 0;
        int startSpace = 0;
        MPI_Recv(&localSize, 1, MPI_INT, 0, TAG_SIZE, MPI_COMM_WORLD, &status);
        MPI_Recv(&startSpace, 1, MPI_INT, 0, TAG_START, MPI_COMM_WORLD, &status);
        double ** res = (double **) malloc(localSize * sizeof(double *));
        for (i = 0; i < localSize; ++i) {
            res[i] = (double *) malloc(JSIZE * sizeof(double));
        }
        for (i = 0; i < localSize; ++i) {
            for (j = 0; j < JSIZE; ++j) {
                res[i][j] = sin(0.00001 * (10 * (i + startSpace) + j));
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

// responsible for printing a into ff
void print_to_file(FILE *ff, double **a) {
    int i, j;
    for (i = 0; i < ISIZE; i++) {
        for (j = 0; j < JSIZE; j++) {
            fprintf(ff, "%f ", a[i][j]);
        }
        fprintf(ff, "\n");
    }
}
