//
// Created by zigal0 on 22.04.2021.
//

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

double FunctionSpace0(double x);

double FunctionTime0(double t);

double FunctionF(double x, double t);

int main(int argc, char **argv) {
    double time_start, time_finish; // measure time
    double time = 1, space = 1; // boundaries
    int t, s; // iterations
    int rowT = 3000, colS = 2520; // quantity of steps (net)
    int localSize; // quantity of iterations in every process
    int dest, src; // destination and source addresses for send & receive
    MPI_Status status; // status for checking delivery
    MPI_Request request; // status of delivery
    int rank = 0, size = 0; // rank of process and quantity of processes
    double **res = NULL; // result

    if (argc == 2) {
        if (1 != sscanf(argv[1], "%d", &rowT)) {
            printf("Error \n");
            return 0;
        }
    }
    if (argc == 3) {
        if (1 != sscanf(argv[1], "%d", &rowT)) {
            printf("Error \n");
            return 0;
        }
    }
    double timeStep = time / rowT, spaceStep = space / colS; // steps

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int startSpace = rank * (colS / size); // start x for iteration in every process
    localSize = colS / size + 1;

    dest = rank + 1;
    if (dest == size) {
        dest = MPI_PROC_NULL;
    }
    src = rank - 1;
    if (src == -1) {
        src = MPI_PROC_NULL;
    }

    if (rank == 0) {
        time_start = MPI_Wtime(); // time of start
        res = (double **) malloc((colS + 1) * sizeof(double *));
        for (s = 0; s <= colS; s++) {
            res[s] = (double *) malloc((rowT + 1) * sizeof(double));
        }
        // left boundary
        for (t = 0; t <= rowT; t++) {
            res[0][t] = FunctionTime0(t * timeStep);
        }
        // bottom boundary
        for (s = 0; s <= colS; s++) {
            res[s][0] = FunctionSpace0(s * spaceStep);
        }
        for (t = 0; t <= rowT - 1; t++) {
            for (s = startSpace + 1; s < startSpace + localSize; s++) {
                res[s][t + 1] = res[s][t] + timeStep * (FunctionF(s * spaceStep, t * timeStep)
                                                        - (res[s][t] - res[s - 1][t]) / spaceStep);
            }
            MPI_Isend(&res[localSize - 1][t + 1], 1, MPI_DOUBLE, dest, 5, MPI_COMM_WORLD, &request);
        }
        // collect all parts of task
        for (int i = 1; i < size; i++) {
            for (int j = (localSize - 1) * i + 1; j < (localSize - 1) * (i + 1) + 1; j++) {
                MPI_Recv(res[j], (rowT + 1), MPI_DOUBLE, i, 6, MPI_COMM_WORLD, &status);
            }
        }
        time_finish = MPI_Wtime(); // time of finish
//        printf("Quantity of processes - %d - time - %f\n",size, time_finish - time_start);
        printf("%f\n", time_finish - time_start);
    } else {
        res = (double **) malloc(localSize * sizeof(double *));
        for (s = 0; s < localSize; s++) {
            res[s] = (double *) malloc((rowT + 1) * sizeof(double));
        }
        int i;
        for (s = 0, i = startSpace; s < localSize; s++, i++) {
            res[s][0] = FunctionSpace0(i * spaceStep);
        }
        for (t = 0; t <= rowT - 1; t++) {
            for (i = startSpace + 1, s = 1; s < localSize; s++, i++) {
                res[s][t + 1] = res[s][t] + timeStep * (FunctionF(i * spaceStep, t * timeStep)
                                                        - (res[s][t] - res[s - 1][t]) / spaceStep);
            }
            // receive and send boundary value
            MPI_Recv(&res[0][t + 1], 1, MPI_DOUBLE, src, 5, MPI_COMM_WORLD, &status);
            MPI_Isend(&res[localSize - 1][t + 1], 1, MPI_DOUBLE, dest, 5, MPI_COMM_WORLD, &request);
        }
        // send part of task to master
        for (s = 1; s < localSize; s++) {
            MPI_Send(res[s], rowT + 1, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD);
        }
    }
    //// result out
//    if (rank == 0) {
//        for (t = rowT; t >= 0; t--) {
//            for (s = 0; s <= colS; s++) {
//                printf("%f\t", res[s][t]);
//            }
//            printf("\n");
//        }
//    }
    // free all structures
    if (rank == 0) {
        for (s = 0; s <= colS; s++) {
            free(res[s]);
        }
    } else {
        for (s = 0; s < localSize; s++) {
            free(res[s]);
        }
    }
    free(res);
    MPI_Finalize();
    return 0;
}

double FunctionSpace0(double x) {
    return 3 * x;
}

double FunctionTime0(double t) {
    return 2 * t;
}

double FunctionF(double x, double t) {
    return x + t;
}