cmake_minimum_required(VERSION 3.0)
project(IntegralOpenMP)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lm")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -lm")
set(CMAKE_CXX_STANDARD 14)

add_definitions(-DBENCHMARK -DSI_PREPARE_COEFFICIENTS)

add_executable(IntegralOpenMP main.cpp functions/arcsin.cpp functions/arcsin.h types.h)
