cmake_minimum_required(VERSION 3.17)
project(hlr.exercise C)

#set(CMAKE_C_STANDARD 90)

add_executable(map hlr.exercise.01/map/map.c)
add_executable(pointer hlr.exercise.01/pointer/pointer.c)

add_executable(simple hlr.exercise.02/simple/simple.c)
add_executable(broken-pde hlr.exercise.02/broken-pde/partdiff-seq.c)
