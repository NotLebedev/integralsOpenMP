cmake_minimum_required(VERSION 3.0)
project(IntegralOpenMP)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lm")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -lm")
set(CMAKE_CXX_STANDARD 14)

add_definitions(-DEXP_PREPARE_COEFFICIENTS -DHS_PREPARE_COEFFICIENTS -DSI_PREPARE_COEFFICIENTS)

add_executable(IntegralOpenMP main.cpp functions/arcsin.cpp functions/arcsin.h types.h functions/heaviside_step.cpp functions/heaviside_step.h functions/exp.cpp functions/exp.h)
