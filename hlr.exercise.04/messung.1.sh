#!/bin/bash

messung=messung.1
mkdir -p $messung;

iterations=600

for iter in `seq 0 2`
do
    for inx in `seq 1 12`
    do
        echo "Messung mit '$inx' threads"
        OMP_NUM_THREADS=$inx ./partdiff-openmp 1 2 512 2 2 $iterations > $messung/$messung.$iter.threads-$inx.txt
    done
done
