//
// Created by zigal0 on 22.04.2021.
//

#include <iostream>
#include <cmath>
#include "mpi.h"
#include <fstream>
#include <filesystem>
#include <cstring>

namespace fs = std::filesystem;

double FunctionSpace0(double x);

double FunctionTime0(double t);

double FunctionF(double x, double t);

int main(int argc, char **argv) {
    std::string path = fs::current_path();
    double a = 2; // coefficient
    double time_start, time_finish; // measure time
    double time = 1, space = 1; // boundaries
    int t, s; // iterations
    int rowT = 1000, colS = 1000; // quantity of steps (net)
    int dest, src; // destination and source addresses for send & receive
    MPI_Status status; // status for checking delivery
    MPI_Request request; // status of delivery
    int rank = 0, size = 0; // rank of process and quantity of processes
    bool sol = false;
    double **res; // result
    // get rowT from terminal
    if (argc > 1) {
        std::stringstream ss(argv[1]);
        if ((ss >> rowT).fail()) {
            // it's not convertible to int
            std::cout << "Error" << std::endl;
            return 0;
        }
        if (argc > 2) {
            if (strcmp(argv[2], "+s") == 0) {
                sol = true;
            }
        }
    }
    double timeStep = time / rowT, spaceStep = space / colS; // steps
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    dest = rank + 1;
    if (dest == size) {
        dest = MPI_PROC_NULL;
    }
    src = rank - 1;
    if (src == -1) {
        src = MPI_PROC_NULL;
    }

    if (rank == 0) {
        int *localSize;  // for uniform distribution between processes
        localSize = (int *) malloc(size * sizeof(int));
        int *startSpace;  // start x for iteration in every process
        startSpace = (int *) malloc(size * sizeof(int));
        int curNum = colS;
        startSpace[0] = 0;
        // quantity of iterations in every process and start
        for (int p = 0; p < size; ++p) {
            if (curNum % (size - p) != 0) {
                localSize[p] = curNum / size + 2;
            } else {
                localSize[p] = curNum / (size - p) + 1;
            }
            curNum -= localSize[p] - 1;
            if (p != 0) {
                MPI_Isend(&localSize[p], 1, MPI_INT, p, 1, MPI_COMM_WORLD, &request);
                MPI_Isend(&startSpace[p], 1, MPI_INT, p, 1, MPI_COMM_WORLD, &request);
            }
            if (p != size - 1) {
                startSpace[p + 1] = startSpace[p] + localSize[p] - 1;
            }
        }
        time_start = MPI_Wtime(); // time of start
        res = (double **) malloc((colS + 1) * sizeof(double *));
        for (s = 0; s <= colS; ++s) {
            res[s] = (double *) malloc((rowT + 1) * sizeof(double));
        }
        // left boundary
        for (t = 0; t <= rowT; ++t) {
            res[0][t] = FunctionTime0(t * timeStep);
        }
        // bottom boundary
        for (s = 0; s < localSize[0]; ++s) {
            res[s][0] = FunctionSpace0(s * spaceStep);
        }
        for (t = 0; t < rowT; ++t) {
            for (s = 1; s < localSize[0]; ++s) {
                res[s][t + 1] = ((2 * spaceStep * timeStep) / (spaceStep + timeStep * a)) *
                                (FunctionF((s + 0.5) * spaceStep, (t + 0.5) * timeStep) +
                                 a * (res[s - 1][t + 1] - res[s][t] + res[s - 1][t]) / (2 * spaceStep) +
                                 (res[s - 1][t] - res[s - 1][t + 1] + res[s][t]) / (2 * timeStep));
            }
            MPI_Isend(&res[localSize[0] - 1][t + 1], 1, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD, &request);
        }
        // collect all parts of task
        for (int i = 1; i < size; ++i) {
            for (int j = startSpace[i] + 1; j < startSpace[i] + localSize[i]; ++j) {
                MPI_Recv(res[j], (rowT + 1), MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status);
            }
        }
        time_finish = MPI_Wtime(); // time of finish
        std::cout << time_finish - time_start << std::endl;
        // result out (csv)
        if (sol) {
            std::ofstream resultFile;
            resultFile.open(path + '/' + "OutPut/solution" + std::to_string(size) + ".csv");
            for (t = 0; t <= rowT; ++t) {
                for (s = 0; s <= colS; ++s) {
                    resultFile << res[s][t];
                    if (s != colS) {
                        resultFile << ',';
                    }
                }
                resultFile << "\n";
            }
            resultFile.close();
        }
        // free all structures except res
        free(startSpace);
        free(localSize);
        for (s = 0; s <= colS; ++s) {
            free(res[s]);
        }

    } else {
        int localSize = 0;
        int startSpace = 0;
        MPI_Recv(&localSize, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&startSpace, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        res = (double **) malloc(localSize * sizeof(double *));
        for (s = 0; s < localSize; ++s) {
            res[s] = (double *) malloc((rowT + 1) * sizeof(double));
        }
        int i;
        // bottom boundary
        for (s = 0, i = startSpace; s < localSize; ++s, ++i) {
            res[s][0] = FunctionSpace0(i * spaceStep);
        }
        for (t = 0; t < rowT; ++t) {
            // receive boundary value
            MPI_Recv(&res[0][t + 1], 1, MPI_DOUBLE, src, 2, MPI_COMM_WORLD, &status);
            for (i = startSpace + 1, s = 1; s < localSize; ++s, ++i) {
                res[s][t + 1] = ((2 * spaceStep * timeStep) / (spaceStep + timeStep * a)) *
                                (FunctionF((i + 0.5) * spaceStep, (t + 0.5) * timeStep) +
                                 a * (res[s - 1][t + 1] - res[s][t] + res[s - 1][t]) / (2 * spaceStep) +
                                 (res[s - 1][t] - res[s - 1][t + 1] + res[s][t]) / (2 * timeStep));
            }
            // send boundary value
            MPI_Isend(&res[localSize - 1][t + 1], 1, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD, &request);
        }
        // send part of task to master
        for (s = 1; s < localSize; ++s) {
            MPI_Send(res[s], rowT + 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
        }
        // free all structures except res
        for (s = 0; s < localSize; ++s) {
            free(res[s]);
        }
    }

    free(res);
    MPI_Finalize();
    return 0;
}

double FunctionSpace0(double x) {
    return cos(M_PI * x);
}

double FunctionTime0(double t) {
    return exp(-t);
}

double FunctionF(double x, double t) {
    return x + t;
}