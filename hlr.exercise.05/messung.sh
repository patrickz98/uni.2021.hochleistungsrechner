#!/bin/bash

# ./partdiff-posix 12 2 512 2 2 10240 > zzz.posix.10240.txt;
OMP_NUM_THREADS=12 ./partdiff-openmp 12 2 512 2 2 10240 > zzz.openmp.10240.txt;
