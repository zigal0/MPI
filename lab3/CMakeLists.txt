cmake_minimum_required(VERSION 3.20)
project(lab3 C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_COMPILER mpicc)

include_directories(/usr/include/x86_64-linux-gnu/mpich)

add_executable(lab3
#        src/consistentMain.c
#        src/parallelMain.c
#        src/parallelMainOther.c
        src/consistent2d.c
#        src/parallel2d.c

        )

target_link_libraries(lab3 ${MPI_C_LIBRARIES} m)
