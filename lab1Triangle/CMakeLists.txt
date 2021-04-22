cmake_minimum_required(VERSION 3.19)
project(lab1Triangle C)

set(CMAKE_C_STANDARD 11)

include_directories(/usr/include/x86_64-linux-gnu/mpich)

set(CMAKE_C_COMPILER "/usr/bin/mpicc")

add_executable(lab1Triangle
#        consistent.c
        parallel.c
        )