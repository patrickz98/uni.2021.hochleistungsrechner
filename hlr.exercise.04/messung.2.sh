#!/bin/bash

name=messung.2
mkdir -p $name;

steps=11
stepSize=$((1024 / $steps))

for iter in `seq 0 2`
do
    for inx in `seq 1 $steps`
    do
        interlines=$(($stepSize * $inx));
        echo "$inx => $interlines";
        OMP_NUM_THREADS=12 ./partdiff-openmp 1 2 $interlines 2 2 100 > $name/$name.$iter.interlines-$interlines.txt
    done
done
