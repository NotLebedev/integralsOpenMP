cmake_minimum_required(VERSION 3.0)
project(IntegralOpenMP)

include_directories(/usr/include/mpi/)
set(CMAKE_C_COMPILER /usr/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/bin/mpic++)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lm")
set(CMAKE_CXX_STANDARD 14)

add_definitions(-DEXP_PREPARE_COEFFICIENTS -DHS_PREPARE_COEFFICIENTS -DSI_PREPARE_COEFFICIENTS)

add_executable(IntegralOpenMP main.cpp functions/exp.cpp functions/heaviside_step.cpp functions/arcsin.cpp messaging.cpp messaging.h)
