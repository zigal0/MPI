//
// Created by zigal0 on 22.04.2021.
//

#include <iostream>
#include <cmath>
#include "mpi.h"
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

double FunctionSpace0(double x);

double FunctionTime0(double t);

double FunctionF(double x, double t);

int main(int argc, char **argv) {
    std::string path = fs::current_path();
    MPI_Init(&argc, &argv);
    double a = 2;
    double time_start, time_finish; // measure time
    double time = 1, space = 1; // boundaries
    int t, s; // iterations
    int rowT = 100, colS = 100; // quantity of steps (net)
    double timeStep = time / rowT, spaceStep = space / colS; // steps
    time_start = MPI_Wtime();

    // result
    double **res;
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
    // Rectangle Scheme
    for (t = 0; t < rowT; t++) {
        for (s = 1; s <= colS; s++) {
            res[s][t + 1] = ((2 * spaceStep * timeStep) / (spaceStep + timeStep * a)) *
                            (FunctionF((s + 0.5) * spaceStep, (t + 0.5) * timeStep) + a *
                                                                                      (res[s - 1][t + 1] - res[s][t] +
                                                                                       res[s - 1][t]) /
                                                                                      (2 * spaceStep) +
                             (res[s - 1][t] - res[s - 1][t + 1] + res[s][t]) / (2 * timeStep));
        }
    }
    // left Angle scheme
//    for (t = 0; t <= rowT - 1; t++) {
//        for (s = 1; s <= colS; s++) {
//            res[s][t + 1] = res[s][t] + timeStep * (FunctionF(s * spaceStep, t * timeStep)
//                                                    - (res[s][t] - res[s - 1][t]) / spaceStep);
//        }
//    }

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
    time_finish = MPI_Wtime();
    // time
    std::cout << time_finish - time_start << std::endl;

    // result out (terminal)
//    for (t = 0; t <= rowT; t++) {
//        for (s = colS; s >= 0; s--) {
//            std::cout << res[s][t] << '\t';
//        }
//        std::cout << std::endl;
//    }
    // result out (csv)
    std::ofstream resultFile;
    resultFile.open(path + '/' + "OutPut/solution(single).csv");
    for (t = 0; t <= rowT; t++) {
        for (s = 0; s <= colS; s++) {
            resultFile << res[s][t];
            if (s != colS) {
                resultFile << ',';
            }
        }
        resultFile << "\n";
    }
    resultFile.close();

    // free all structures
    for (s = 0; s <= colS; s++) {
        free(res[s]);
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
