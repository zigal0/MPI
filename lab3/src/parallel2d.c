//
// Created by zigal0 on 17.11.2021.
//
// There are 2 realization:
// First realization:
// 1) parallel outer loop
// 2) parallel inner loop
//
// Second realization:
// 1) parallel outer loop
// 2) parallel outer loop
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

// size of array to calculate
#define ISIZE 1000
#define JSIZE 1000

// tags for send/receive functions
#define TAG_SIZE 1
#define TAG_POS 2
#define TAG_RESULT 3
#define TAG_CONTINUE 4
#define TAG_FINAL 5
#define TAG_ARRAY 6

void compute_solo(double **a);

void print_ff(FILE *ff, double **a);

void distribute_1st_cycle_master(double **a, int size);

void compute_1st_cycle_slave();

void distribute_2nd_cycle_master(double **a, int size);

void compute_2nd_cycle_slave();

void distribute_2nd_cycle_master2(double **a, int size);

void compute_2nd_cycle_slave2();

int main(int argc, char **argv) {
    // variables
    int rank = 0, size = 0;         // rank of process and quantity of processes
    int i;                          // for iteration
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // check quantity of processes
    if (rank == 0) {
        if (size < 2) {
            fprintf(stderr, "Error: number of processes must be two or more.\n");
            exit(-1);
        }
    }

    // master
    if (rank == 0) {
        // variables
        FILE *ff;                                            // file to write
        double time_start, time_finish;                      // for measuring time
        double **a = malloc(ISIZE * sizeof(double *));  // array to calculate
        for (i = 0; i < ISIZE; i++) {
            a[i] = malloc(JSIZE * sizeof(double));
        }
        // parallel
        // method # 1
        time_start = MPI_Wtime();
        distribute_1st_cycle_master(a, size);
        distribute_2nd_cycle_master(a, size);
        time_finish = MPI_Wtime();
        printf("Method # 1: execution time -  %0.6f\n", time_finish - time_start);

        // print to txt
        ff = fopen("output/resultP2dM1.txt", "w");
        print_ff(ff, a);
        fclose(ff);

        // method # 2
//         check quantity of processes
        if (size > 9) {
            fprintf(stderr, "Impossible to execute method # 2 because number of processes > 9.\n");
        } else {
            time_start = MPI_Wtime();
            distribute_1st_cycle_master(a, size);
            distribute_2nd_cycle_master2(a, size);
            time_finish = MPI_Wtime();
            printf("Method # 2: execution time -  %0.6f\n", time_finish - time_start);

            // print to txt
            ff = fopen("output/resultP2dM2.txt", "w");
            print_ff(ff, a);
            fclose(ff);
        }

        // consistent
        time_start = MPI_Wtime();
        compute_solo(a);
        time_finish = MPI_Wtime();
        printf("Consistent: execution time - %0.6f\n", time_finish - time_start);

        // print to txt
        ff = fopen("output/resultS2d.txt", "w");
        print_ff(ff, a);
        fclose(ff);

        // free all memory
        for (i = 0; i < ISIZE; i++) {
            free(a[i]);
        }
        free(a);

        // slaves
    } else {
        // method # 1
        compute_1st_cycle_slave();
        compute_2nd_cycle_slave();

        // method # 2
        compute_1st_cycle_slave();
        compute_2nd_cycle_slave2();
    }

    // waiting for all processes
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

void compute_solo(double **a) {
    int i, j;
    for (i = 0; i < ISIZE; i++) {
        for (j = 0; j < JSIZE; j++) {
            a[i][j] = 10 * i + j;
        }
    }
    for (i = 8; i < ISIZE; i++) {
        for (j = 0; j < JSIZE - 3; j++) {
            a[i][j] = sin(0.00001 * a[i - 8][j + 3]);
        }
    }
}

void print_ff(FILE *ff, double **a) {
    for (int i = 0; i < ISIZE; i++) {
        for (int j = 0; j < JSIZE; j++) {
            fprintf(ff, "%f ", a[i][j]);
        }
    }
}

void distribute_1st_cycle_master(double **a, int size) {
    // variables
    MPI_Status status;                                          // status for checking delivery
    MPI_Request request;                                        // status of delivery
    int numSlaves = size - 1;                                   // number of slaves
    int curNum = ISIZE;                                         // current number of columns to distribute
    int p;                                                      // for iteration
    int *localSize = (int *) malloc(size * sizeof(int));   // for uniform distribution between processes
    int *localPos = (int *) malloc(size * sizeof(int));    // start position for iteration in every process
    localPos[0] = 0;

    // calculate localPos & localSize and send tasks
    for (p = 0; p < numSlaves; ++p) {
        if (curNum % (numSlaves - p) != 0) {
            localSize[p] = curNum / (numSlaves - p) + 1;
        } else {
            localSize[p] = curNum / (numSlaves - p);
        }
        curNum -= localSize[p];
        MPI_Isend(&localSize[p], 1, MPI_INT, p + 1, TAG_SIZE, MPI_COMM_WORLD, &request);
        MPI_Isend(&localPos[p], 1, MPI_INT, p + 1, TAG_POS, MPI_COMM_WORLD, &request);
        if (p != numSlaves - 1) {
            localPos[p + 1] = localPos[p] + localSize[p];
        }
    }

    // collect all parts of task
    int curPos = 0;
    for (p = 0; p < numSlaves; ++p) {
        for (int j = localPos[p]; j < localPos[p] + localSize[p]; ++j) {
            MPI_Recv(a[curPos], JSIZE, MPI_DOUBLE, p + 1, TAG_RESULT, MPI_COMM_WORLD, &status);
            curPos++;
        }
    }

    // free all memory
    free(localSize);
    free(localPos);
}

void compute_1st_cycle_slave() {
    // variables
    MPI_Status status;          // status for checking delivery
    int i, j;                   // for iteration
    int localSize = 0;          // size of batch
    int localPos = 0;           // start position in general array

    // get required parameters from master
    MPI_Recv(&localSize, 1, MPI_INT, 0, TAG_SIZE, MPI_COMM_WORLD, &status);
    MPI_Recv(&localPos, 1, MPI_INT, 0, TAG_POS, MPI_COMM_WORLD, &status);

    // memory for result
    double **res = (double **) malloc(localSize * sizeof(double *));
    for (i = 0; i < localSize; ++i) {
        res[i] = (double *) malloc(JSIZE * sizeof(double));
    }

    // calculate result
    for (i = 0; i < localSize; ++i) {
        for (j = 0; j < JSIZE; ++j) {
            res[i][j] = 10 * (i + localPos) + j;
        }
    }

    // send part of task to master
    for (i = 0; i < localSize; ++i) {
        MPI_Send(res[i], JSIZE, MPI_DOUBLE, 0, TAG_RESULT, MPI_COMM_WORLD);
    }

    // free all memory
    for (i = 0; i < localSize; ++i) {
        free(res[i]);
    }
    free(res);
}

void distribute_2nd_cycle_master(double **a, int size) {
    // variables
    MPI_Request request;                                        // status of delivery
    MPI_Status status;                                          // status for checking delivery
    int i, p;                                                   // for iteration
    int numSlaves = size - 1;                                          // number of slaves
    int curNum = JSIZE - 3;                                     // current number of rows to distribute
    int *localSize = (int *) malloc(numSlaves * sizeof(int));      // for uniform distribution between processes
    int *localPos = (int *) malloc(numSlaves * sizeof(int));       // start pos for iteration in every process
    localPos[0] = 0;

    // calculate localPos & localSize
    for (p = 0; p < numSlaves; ++p) {
        if (curNum % (numSlaves - p) != 0) {
            localSize[p] = curNum / (numSlaves - p) + 1;
        } else {
            localSize[p] = curNum / (numSlaves - p);
        }
        curNum -= localSize[p];
        if (p != numSlaves - 1) {
            localPos[p + 1] = localPos[p] + localSize[p];
        }
    }

    // send localSizes to the processes
    for (p = 0; p < numSlaves; ++p) {
        MPI_Isend(&localSize[p], 1, MPI_INT, p + 1, TAG_SIZE, MPI_COMM_WORLD, &request);
    }

    // send & receive tasks
    for (i = 8; i < ISIZE; ++i) {
        for (p = 0; p < numSlaves; ++p) {
            if (i == ISIZE - 1) {
                MPI_Isend(a[i - 8] + localPos[p] + 3, localSize[p], MPI_DOUBLE, p + 1, TAG_FINAL,
                          MPI_COMM_WORLD, &request);
            } else {
                MPI_Isend(a[i - 8] + localPos[p] + 3, localSize[p], MPI_DOUBLE, p + 1, TAG_CONTINUE,
                          MPI_COMM_WORLD, &request);
            }
        }
        for (p = 0; p < numSlaves; ++p) {
            MPI_Recv(a[i] + localPos[p], localSize[p], MPI_DOUBLE, p + 1, TAG_RESULT, MPI_COMM_WORLD,
                     &status);
        }
    }

    // free all memory
    free(localSize);
    free(localPos);
}

void compute_2nd_cycle_slave() {
    // variables
    int localSize = 0;          // size of batch
    MPI_Status status;          // status for checking delivery

    // get size of batch
    MPI_Recv(&localSize, 1, MPI_INT, 0, TAG_SIZE, MPI_COMM_WORLD, &status);

    // memory for result
    double *res = malloc(localSize * sizeof(double));

    // calculate result
    while (1) {
        MPI_Recv(res, localSize, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        for (int i = 0; i < localSize; ++i) {
            res[i] = sin(0.00001 * res[i]);
        }
        MPI_Send(res, localSize, MPI_DOUBLE, 0, TAG_RESULT, MPI_COMM_WORLD);
        if (status.MPI_TAG == TAG_FINAL) {
            break;
        }
    }
    // free all memory
    free(res);
}

void distribute_2nd_cycle_master2(double **a, int size) {
    // variables
    int i;                          // for iteration
    int pos = 0;                    // current position
    int numSlaves = size - 1;       // number of slaves
    int colsParalleled = 0;         // current number of paralleled columns
    int numProcRunning;             // number of running processes
    int numProcFinished;            // number of finished processes
    int restriction = 8;            // max number of parallel processes (in our case 8)
    MPI_Status status;              // status for checking delivery

    // distribution
    while (colsParalleled < restriction) {
        numProcRunning = 0;
        numProcFinished = 0;
        for (i = 0; i < numSlaves; i++) {
            if (colsParalleled < restriction) {
                pos = restriction + colsParalleled;
                MPI_Send(&pos, 1, MPI_INT, i + 1, TAG_POS, MPI_COMM_WORLD);
                MPI_Send(a[pos], JSIZE, MPI_DOUBLE, i + 1, TAG_ARRAY, MPI_COMM_WORLD);
                MPI_Send(a[pos - restriction], JSIZE, MPI_DOUBLE, i + 1, TAG_ARRAY, MPI_COMM_WORLD);
                numProcRunning++;
                colsParalleled++;
            }
        }
        while (numProcFinished < numProcRunning) {
            MPI_Recv(&pos, 1, MPI_INT, MPI_ANY_SOURCE, TAG_POS, MPI_COMM_WORLD, &status);
            MPI_Recv(a[pos], JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY, MPI_COMM_WORLD, &status);
            if (pos < ISIZE - restriction) {
                pos += restriction;
                MPI_Send(&pos, 1, MPI_INT, status.MPI_SOURCE, TAG_POS, MPI_COMM_WORLD);
                MPI_Send(a[pos], JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY, MPI_COMM_WORLD);
            } else {
                numProcFinished++;
            }
        }
        if ((ISIZE - restriction) / restriction == 0) {
            break;
        }
    }

    // stop all slaves
    for (i = 0; i < numSlaves; i++) {
        MPI_Send(&pos, 1, MPI_INT, i + 1, TAG_FINAL, MPI_COMM_WORLD);
    }
}

void compute_2nd_cycle_slave2() {
    // variables
    double *a_i = malloc(JSIZE * sizeof(double));       // current column to calculate
    double *a_i_8 = malloc(JSIZE * sizeof(double));     // prev column
    int pos = 0;                                            // current position
    int j;                                                  // for iteration
    double *save;                                           // for exchange
    MPI_Status status;                                      // status for checking delivery

    while (1) {

        // receive position (check for break)
        MPI_Recv(&pos, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == TAG_FINAL) {
            break;
        }
        // receive new column & receive addition column if 1st time
        MPI_Recv(a_i, JSIZE, MPI_DOUBLE, 0, TAG_ARRAY, MPI_COMM_WORLD, &status);
        if (pos < 16) {
            MPI_Recv(a_i_8, JSIZE, MPI_DOUBLE, 0, TAG_ARRAY, MPI_COMM_WORLD, &status);
        }
        // calculation
        for (j = 0; j < JSIZE - 3; j++) {
            a_i[j] = sin(0.00001 * a_i_8[j + 3]);
        }
        // send current pos and calculated column
        MPI_Send(&pos, 1, MPI_INT, 0, TAG_POS, MPI_COMM_WORLD);
        MPI_Send(a_i, JSIZE, MPI_DOUBLE, 0, TAG_ARRAY, MPI_COMM_WORLD);
        save = a_i_8;
        a_i_8 = a_i;
        a_i = save;
    }

    // free all memory
    free(a_i);
    free(a_i_8);
}