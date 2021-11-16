//
// Created by zigal0 on 25.04.2021.
//

#include <pthread.h>
#include <iostream>
#include <cmath>
#include <chrono>

typedef struct thread_data_t {
    double start, finish;
} thread_data_t;

pthread_mutex_t lock;
double res = 0;
int power = -15;
int steps;

double func(double x);

void *mythread(void *arg);

int main() {
    double leftB = 0.001, rightB = 1; // boundaries
    // data
//    std::cout << "Integral of sin(1/x) in the range x = leftB..rightB with accuracy (e = 10 ^(p)" << std::endl;
//    std::cout << "The interval is divided by the intervals between zeros with the same numbers of steps: " << std::endl;
//    std::cout << "Enter range :" << std::endl;
//    std::cout << "leftB = ";
//    std::cin >> leftB;
//    std::cout << "rightB = ";
//    std::cin >> rightB;
//    std::cout << "Enter the power for desired accuracy (e = 10 ^(p)), p = ";
//    std::cin >> power;

    int num = (int) (1 / (M_PI * leftB));
    double cur = leftB; // current x
    int i, rc;
    pthread_t thr[num + 1];
    thread_data_t thr_data[num + 1];

    std::cout.precision(-power); // a number of symbols after comma
    pthread_mutex_init(&lock, nullptr); // initiate mutex

    auto begin = std::chrono::steady_clock::now(); // start time
    steps = (int) ((rightB - 1 / M_PI) / pow((6480 * pow(1 / M_PI, 8) * pow(10, power)) / (rightB - 1 / M_PI), 0.25));
    for (i = num; i > 0; --i) {
        thr_data[i].start = cur;
        cur = 1 / (M_PI * i);
        thr_data[i].finish = 1 / (M_PI * i);
        if ((rc = pthread_create(&thr[i], nullptr, mythread, &thr_data[i]))) {
            fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
            return EXIT_FAILURE;
        }
    }
    thr_data[0].start = cur;
    thr_data[0].finish = rightB;
    if ((rc = pthread_create(&thr[0], nullptr, mythread, &thr_data[i]))) {
        fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
        return EXIT_FAILURE;
    }
    // waiting for all threads
    for (i = 0; i <= num; ++i) {
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
    double step = (data->finish - data->start) / steps;
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