//
// Created by zigal0 on 22.04.2021.
//

#include <stdio.h>
#include <malloc.h>
#include "mpi.h"

double FunctionSpace0(double x);

double FunctionTime0(double t);

double FunctionF(double x, double t);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    double time_start, time_finish; // measure time
    double time = 1, space = 1; // boundaries
    int t, s; // iterations
    int rowT = 3000, colS = 2520; // quantity of steps (net)
    double timeStep = time / rowT, spaceStep = space / colS; // steps
    time_start = MPI_Wtime();
    // result
    double **res = NULL;
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
    // left Angle scheme
    for (t = 0; t <= rowT - 1; t++) {
        for (s = 1; s <= colS; s++) {
            res[s][t + 1] = res[s][t] + timeStep * (FunctionF(s * spaceStep, t * timeStep)
                                                    - (res[s][t] - res[s - 1][t]) / spaceStep);
        }
    }
    // time
    time_finish = MPI_Wtime();
    printf("%f\n", time_finish - time_start);
    // cross scheme + left Angle 1-st step and right boundary
//    for (s = 1; s <= colS; s++) {
//        res[s][1] = res[s][0] + timeStep * (FunctionF(s * spaceStep, 0)
//                - (res[s][0] - res[s - 1][0]) / spaceStep);
//    }
//    for (t = 1; t < rowT; t++) {
//        for (s = 1; s < colS; s++) {
//            res[s][t + 1] = res[s][t - 1] + 2 * timeStep * (FunctionF(s * spaceStep, t * timeStep)
//                    - (res[s + 1][t] - res[s - 1][t]) / (2 * spaceStep));
//        }
//        res[s][t + 1] = res[s][t] + timeStep * (FunctionF(s * spaceStep, t * timeStep)
//                                            - (res[s][t] - res[s - 1][t]) / spaceStep);
//    }
    // result out
//    for (t = rowT; t >= 0; t--) {
//        for (s = 0; s <= colS; s++) {
//            printf("%f\t", res[s][t]);
//        }
//        printf("\n");
//    }


    // free all structures
    for (s = 0; s <= colS; s++) {
        free(res[s]);
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