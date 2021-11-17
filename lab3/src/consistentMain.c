//
// Created by zigal0 on 16.11.2021.
//
// Consistent realization of main task with measuring time
//

#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <malloc.h>

#define ISIZE 1000
#define JSIZE 1000

int main(int argc, char **argv) {
    int i, j;
    double **a = malloc(ISIZE * sizeof(double *));
    for (i = 0; i < ISIZE; i++) {
        a[i] = malloc(JSIZE * sizeof(double));
    }
    FILE *ff;
    double time_start, time_finish; // measure time
    MPI_Init(&argc, &argv);
    time_start = MPI_Wtime();
    for (i = 0; i < ISIZE; i++) {
        for (j = 0; j < JSIZE; j++) {
            a[i][j] = 10 * i + j;
        }
    }
    for (i = 0; i < ISIZE; i++) {
        for (j = 0; j < JSIZE; j++) {
            a[i][j] = sin(0.00001 * a[i][j]);
        }
    }
    time_finish = MPI_Wtime();
    printf("execution time: %0.6f\n", time_finish - time_start);
    ff = fopen("output/resultCM.txt", "w");
    for (i = 0; i < ISIZE; i++) {
        for (j = 0; j < JSIZE; j++) {
            fprintf(ff, "%f ", a[i][j]);
        }
        fprintf(ff, "\n");
    }
    fclose(ff);
    for (i = 0; i < ISIZE; ++i) {
        free(a[i]);
    }
    free(a);
    return 0;
}
