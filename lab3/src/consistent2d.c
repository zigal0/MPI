//
// Created by zigal0 on 16.11.2021.
//
// Consistent realization of 2d task with measuring time
//

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <mpi.h>

#define ISIZE 1000
#define JSIZE 1000

int main(int argc, char **argv) {
    int i, j;
    double time_start, time_finish;
    double **a = malloc(ISIZE * sizeof(double *));
    for (i = 0; i < ISIZE; i++) {
        a[i] = malloc(JSIZE * sizeof(double));
    }
    MPI_Init(&argc, &argv);
    time_start = MPI_Wtime();
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
    time_finish = MPI_Wtime();
    printf("execution time: %0.6f\n", time_finish - time_start);
    FILE *ff = fopen("output/resultC2d.txt", "w");
    for (i = 0; i < ISIZE; i++) {
        for (j = 0; j < JSIZE; j++) {
            fprintf(ff, "%f ", a[i][j]);
        }
    }
    fclose(ff);
    for (i = 0; i < ISIZE; ++i) {
        free(a[i]);
    }
    free(a);
    return 0;
}

