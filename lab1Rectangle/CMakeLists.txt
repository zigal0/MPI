cmake_minimum_required(VERSION 3.19)
project(lab1Rectangle)

set(CMAKE_CXX_STANDARD 20)

include_directories(/usr/include/x86_64-linux-gnu/mpich)

set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")

add_executable(lab1Rectangle
#        consistent.cpp
        parallel.cpp
        )
