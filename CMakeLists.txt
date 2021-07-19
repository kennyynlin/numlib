cmake_minimum_required(VERSION 3.19)
project(numlib)

set(CMAKE_CXX_STANDARD 11)
find_package(OpenMP)
add_library(numlib SHARED numlib.cpp numlib.h)

if(OpenMP_CXX_FOUND)
    target_link_libraries(numlib PUBLIC OpenMP::OpenMP_CXX)
endif()