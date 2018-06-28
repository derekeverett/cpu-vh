cmake_minimum_required(VERSION 2.8)
project(cpuvh)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

# Require MPI and GSL for this project:
find_package(GSL REQUIRED)

set(CMAKE_CXX_FLAGS "-O3 -fopenmp -march=native -fopt-info-vec")

include_directories(${MPI_INCLUDE_PATH} ${GSL_INCLUDE_DIR})
add_subdirectory (rhic/src)

