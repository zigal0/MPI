//
// Created by zigal0 on 24.04.2021.
//

#include <pthread.h>
#include <iostream>
#include <cmath>
#include <chrono>

#define NUM_THREADS 13

typedef struct thread_data_t {
    double start, finish;
} thread_data_t;

pthread_mutex_t lock;
double res = 0;
int power = -15;

double func(double x);

void *mythread(void *arg);

int main() {
    double leftB = 0.001, rightB = 1; // boundaries
    // data
    std::cout << "Integral of sin(1/x) in the range x = leftB..rightB with accuracy (e = 10 ^(p)" << std::endl;
    std::cout << "The interval is divided by parts like interval / 2^(NUM_THREADS - 1): " << std::endl;
    std::cout << "Enter range :" << std::endl;
    std::cout << "leftB = ";
    std::cin >> leftB;
    std::cout << "rightB = ";
    std::cin >> rightB;
    std::cout << "Enter the power for desired accuracy (e = 10 ^(p)), p = ";
    std::cin >> power;

    double cur = leftB; // current x
    double interval = rightB - leftB;
    double variableStep = interval / pow(2, NUM_THREADS - 1);
    int i, rc;
    pthread_t thr[NUM_THREADS];
    thread_data_t thr_data[NUM_THREADS];

    std::cout.precision(-power); // a number of symbols after comma
    pthread_mutex_init(&lock, nullptr); // initiate mutex

    auto begin = std::chrono::steady_clock::now(); // start time

    for (i = 0; i < NUM_THREADS; ++i) {
        if (i != 0) {
            variableStep = interval / pow(2, NUM_THREADS - i);
        }
        thr_data[i].start = cur;
        cur += variableStep;
        thr_data[i].finish = cur;

        if ((rc = pthread_create(&thr[i], nullptr, mythread, &thr_data[i]))) {
            fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
            return EXIT_FAILURE;
        }
    }
    // waiting for all threads
    for (i = 0; i < NUM_THREADS; ++i) {
        pthread_join(thr[i], nullptr);
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
    auto *data = (thread_data_t *) arg;
    double step = pow((6480 * pow(data->start, 8) * pow(10, power)) / (data->finish - data->start), 0.25);
    int steps = (int) ((data->finish - data->start) / step) + 1;
    step = (data->finish - data->start) / steps;
    double x = data->start;
    for (int i = 0; i < steps; ++i) {
        part += step * (func(x) + 3 * func(x + step / 3) + 3 * func(x + 2 * step / 3) + func(x + step)) / 8;
        x += step;
    }
    pthread_mutex_lock(&lock);
    res += part;
    pthread_mutex_unlock(&lock);
    pthread_exit(nullptr);
}

double func(double x) {
    return sin(1 / x);
}

