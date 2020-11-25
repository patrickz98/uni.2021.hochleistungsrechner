#!/bin/bash

name=messung.2
mkdir -p $name;

steps=11
stepSize=$((1024 / $steps))

for inx in `seq 1 $steps`
do
    interlines=$(($stepSize * $inx));
    echo "$inx => $interlines";
    OMP_NUM_THREADS=12 ./partdiff 1 2 $interlines 2 2 200 > $name/$name.interlines-$interlines.txt
done
