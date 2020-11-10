cmake_minimum_required(VERSION 3.17)
project(hr.exercise.02 C)

set(CMAKE_C_STANDARD 99)

add_executable(simple simple/simple.c)
add_executable(broken-pde broken-pde/partdiff-seq.c)