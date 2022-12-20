cmake_minimum_required(VERSION 3.21)
project(gauss)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(
        bprinter/
)

add_executable(gauss test.cpp)
