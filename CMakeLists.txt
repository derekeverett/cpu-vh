cmake_minimum_required(VERSION 2.8)
#set(CMAKE_CXX_COMPILER "g++")
project(cpu-vh)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

# Require GSL for this project:
find_package(GSL REQUIRED)

# Require libconfig for this project:
#find_package(Config REQUIRED)

#g++ flags
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -std=c++11")

#intel flags
#set(CMAKE_CXX_FLAGS "-O3 -qopenmp -std=c++11")

include_directories(${GSL_INCLUDE_DIR})
#include_directories(${Config_INCLUDE_DIR})

#header files
include_directories(rhic/include)

#freezeout header files
include_directories(freezeout)
include_directories(freezeout/cornelius-c++-1.3)

add_subdirectory (rhic/src)
