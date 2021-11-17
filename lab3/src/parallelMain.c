//
// Created by zigal0 on 16.11.2021.
//
// First realization parallel main task (master distribute only)
// There are 2 addition way: same structure & changed structure
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

// size of array to calculate
#define ISIZE 1000
#define JSIZE 1000

// tags for send/receive functions
#define TAG_STOP 1
#define TAG_PLACE 2
#define TAG_ARRAY 3

void distribute_cycle_master(double **a, int world_size);

void compute_1st_cycle_slave(double *a_i);

void compute_2nd_cycle_slave(double *a_i);

void compute_both_cycles_slave(double *a_i);

void print_to_file(FILE *ff, double **a);

void compute_solo(double **a);

int main(int argc, char **argv) {
    // variables
    int rank = 0, size = 0;         // rank of process and quantity of processes
    int i;                          // iterations

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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
        FILE *ff;                       // file to write
        double time_start, time_finish; // measure time
        double **a = malloc(ISIZE * sizeof(double *));
        for (i = 0; i < ISIZE; i++) {
            a[i] = malloc(JSIZE * sizeof(double));
        }

        // parallel
        // testing parallel structure with 2 cycles
        time_start = MPI_Wtime();
        distribute_cycle_master(a, size);
        distribute_cycle_master(a, size);
        time_finish = MPI_Wtime();
        printf("Parallel: 2 cycles execution time - %0.6f\n", time_finish - time_start);

        // print ro txt
        ff = fopen("output/resultPMM1.txt", "w");
        print_to_file(ff, a);
        fclose(ff);

        // testing parallel structure with 1 cycle
        time_start = MPI_Wtime();
        distribute_cycle_master(a, size);
        time_finish = MPI_Wtime();
        printf("Parallel: 1 cycle execution time - %0.6f\n", time_finish - time_start);

        // print ro txt
        ff = fopen("output/resultPMM2.txt", "w");
        print_to_file(ff, a);
        fclose(ff);

        // consistent
        time_start = MPI_Wtime();
        compute_solo(a);
        time_finish = MPI_Wtime();
        printf("Consistent: execution time - %0.6f\n", time_finish - time_start);

        // print to txt
        ff = fopen("output/resultSM.txt", "w");
        print_to_file(ff, a);
        fclose(ff);

        for (i = 0; i < ISIZE; i++) {
            free(a[i]);
        }
        free(a);

    // slaves
    } else {
        double *a_i = malloc(JSIZE * sizeof(double));
        // for 2 cycles
        compute_1st_cycle_slave(a_i);
        compute_2nd_cycle_slave(a_i);
        // for 1 cycles
        compute_both_cycles_slave(a_i);
        free(a_i);
    }

    // waiting for all processes
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

// main function of master, responsible for assignment distribution and gathering results
void distribute_cycle_master(double **a, int size) {
    // variables
    int pos = 0;                    // current position
    int i;                          // for iteration
    MPI_Status status;              // status for checking delivery

    for (i = 0; i < size - 1; i++) {
        MPI_Send(&i, 1, MPI_INT, i + 1, TAG_PLACE, MPI_COMM_WORLD);
        MPI_Send(a[i], JSIZE, MPI_DOUBLE, i + 1, TAG_ARRAY, MPI_COMM_WORLD);
    }

    for (i = size - 1; i < ISIZE; i++) {
        MPI_Recv(&pos, 1, MPI_INT, MPI_ANY_SOURCE, TAG_PLACE, MPI_COMM_WORLD, &status);
        MPI_Recv(a[pos], JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY, MPI_COMM_WORLD, &status);
        MPI_Send(&i, 1, MPI_INT, status.MPI_SOURCE, TAG_PLACE, MPI_COMM_WORLD);
        MPI_Send(a[i], JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY, MPI_COMM_WORLD);
    }

    for (i = 0; i < size - 1; i++) {
        MPI_Recv(&pos, 1, MPI_INT, MPI_ANY_SOURCE, TAG_PLACE, MPI_COMM_WORLD, &status);
        MPI_Recv(a[pos], JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY, MPI_COMM_WORLD, &status);
        MPI_Send(&pos, 1, MPI_INT, status.MPI_SOURCE, TAG_STOP, MPI_COMM_WORLD);
    }
}

// 1st salve function that calculate 1st cycle of main task
void compute_1st_cycle_slave(double *a_i) {
    // variables
    int pos = 0;            // current position
    int j;                  // for iteration
    MPI_Status status;      // status for checking delivery

    while (1) {
        // receive data
        MPI_Recv(&pos, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == TAG_STOP) {
            break;
        }
        MPI_Recv(a_i, JSIZE, MPI_DOUBLE, 0, TAG_ARRAY, MPI_COMM_WORLD, &status);
        // calculation
        for (j = 0; j < JSIZE; j++) {
            a_i[j] = 10 * pos + j;
        }
        // send data
        MPI_Send(&pos, 1, MPI_INT, 0, TAG_PLACE, MPI_COMM_WORLD);
        MPI_Send(a_i, JSIZE, MPI_DOUBLE, 0, TAG_ARRAY, MPI_COMM_WORLD);
    }
}

// 2nd salve function that calculate 2nd cycle of main task
void compute_2nd_cycle_slave(double *a_i) {
    // variables
    int pos = 0;            // current position
    int j;                  // for iteration
    MPI_Status status;      // status for checking delivery

    while (1) {
        // receive data
        MPI_Recv(&pos, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == TAG_STOP) {
            break;
        }
        MPI_Recv(a_i, JSIZE, MPI_DOUBLE, 0, TAG_ARRAY, MPI_COMM_WORLD, &status);
        // calculation
        for (j = 0; j < JSIZE; j++) {
            a_i[j] = sin(0.00001 * a_i[j]);
        }
        // send data
        MPI_Send(&pos, 1, MPI_INT, 0, TAG_PLACE, MPI_COMM_WORLD);
        MPI_Send(a_i, JSIZE, MPI_DOUBLE, 0, TAG_ARRAY, MPI_COMM_WORLD);
    }
}

// slave function = union of 2 previous functions (excludes additional resend)
void compute_both_cycles_slave(double *a_i) {
    // variables
    int pos = 0;            // current position
    int j;                  // for iteration
    MPI_Status status;      // status for checking delivery

    while (1) {
        // receive data
        MPI_Recv(&pos, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == TAG_STOP) {
            break;
        }
        MPI_Recv(a_i, JSIZE, MPI_DOUBLE, 0, TAG_ARRAY, MPI_COMM_WORLD, &status);
        // calculation
        for (j = 0; j < JSIZE; j++) {
            a_i[j] = sin(0.00001 * (10 * pos + j));
        }
        // send data
        MPI_Send(&pos, 1, MPI_INT, 0, TAG_PLACE, MPI_COMM_WORLD);
        MPI_Send(a_i, JSIZE, MPI_DOUBLE, 0, TAG_ARRAY, MPI_COMM_WORLD);
    }
}

// responsible for printing a into ff
void print_to_file(FILE *ff, double **a) {
    for (int i = 0; i < ISIZE; i++) {
        for (int j = 0; j < JSIZE; j++) {
            fprintf(ff, "%f ", a[i][j]);
        }
        fprintf(ff, "\n");
    }
}

// structure of main task (no parallel here)
void compute_solo(double **a) {
    int i, j;
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
}
