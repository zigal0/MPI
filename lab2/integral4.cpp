//
// Created by zigal0 on 12.05.2021.
//

#include <pthread.h>
#include <iostream>
#include <cmath>
#include <chrono>

#define NUM_THREADS 8
#define INTERVALS 319
#define IT 1000000

typedef struct thread_data_t {
    double start;
    double finish;
} thread_data_t;

pthread_mutex_t lock;
double res = 0;
int power = -12;

double func(double x);

void *mythread(void *arg);

int main() {
    double leftB = 0.001, rightB = 1; // boundaries
    // data
    std::cout << "Integral of sin(1/x) in the range x = 0.001..1 with accuracy (e = 10 ^(p)" << std::endl;
    std::cout << "Every thread has the same quantity of steps" << std::endl;
    std::cout << "Enter the power for desired accuracy (e = 10 ^(p)), p = ";
    std::cin >> power;

    pthread_t thr[NUM_THREADS];
    thread_data_t thr_data[NUM_THREADS];

    std::cout.precision(-power); // a number of symbols after comma
    pthread_mutex_init(&lock, nullptr); // initiate mutex

    auto begin = std::chrono::steady_clock::now(); // start time
    int curIntervals = INTERVALS;
    double curX = rightB;
    int curK;
    int rc;
    for (int i = 0; i < NUM_THREADS; ++i) {
        if (curIntervals % (NUM_THREADS - i) != 0) {
            curK = curIntervals / (NUM_THREADS - i) + 1;
        } else {
            curK = curIntervals / (NUM_THREADS - i);
        }
        curIntervals -= curK;
        thr_data[i].finish = curX;
        thr_data[i].start = 1 / (M_PI * (INTERVALS - curIntervals));
        curX = thr_data[i].start;
        if (i == NUM_THREADS - 1) {
            thr_data[i].start = leftB;
        }

        if ((rc = pthread_create(&thr[i], nullptr, mythread, &thr_data[i]))) {
            fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
            return EXIT_FAILURE;
        }
    }
    // waiting for all threads
    for (unsigned long i : thr) {
        pthread_join(i, nullptr);
    }

    auto end = std::chrono::steady_clock::now(); // finish time
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    // results
    std::cout << "Execution time : " << elapsed_ms.count() << " ns\n";
    std::cout << "Integral of sin(1/x) in the range x = " << leftB << ".." << rightB << " is " << res
              << " with accuracy e = 10^(" << power << ")." << std::endl;
    return EXIT_SUCCESS;
}

void *mythread(void *arg) {
    double part = 0;
    double curRes = 0;
    auto *data = (thread_data_t *) arg;
    double x = data->start;
    double bigStep = (data->finish - data->start) / IT;
    double curStart = x;
    double curFinish = x + bigStep;
    double step = sqrt(12 * pow(10, power) / (bigStep *
                                              (2 * pow(curStart, -3) + pow(curStart, -4))));

    int steps = ceil((curFinish - curStart) / step);
    for (int j = 0; j < IT; j++) {
        step = (curFinish - curStart) / steps;
        for (int i = 0; i < steps; ++i) {
            part += step * (func(x) + func(x + step)) / 2;
            x += step;
        }
        curStart = curFinish;
        curFinish += bigStep;
        curRes += part;
    }
    pthread_mutex_lock(&lock);
    res += part;
    pthread_mutex_unlock(&lock);
    pthread_exit(nullptr);
}

double func(double x) {
    return sin(1 / x);
}

