cmake_minimum_required(VERSION 3.19)
project(lab2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -g")

set(SOURCE_FILES
#        integral1.cpp
#        integral2.cpp
#        integral3.cpp
        integral4.cpp
        )

add_executable(lab2 ${SOURCE_FILES})

find_library(PTHREAD_LIBRARY pthread)

target_link_libraries(lab2 ${PTHREAD_LIBRARY})