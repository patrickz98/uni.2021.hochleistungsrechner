#!/bin/bash

mkdir -p messung.1;

iterations=600

for inx in `seq 1 12`
do
    echo "Messung mit '$inx' threads"
    OMP_NUM_THREADS=$inx ./partdiff 1 2 512 2 2 $iterations > messung.1/messung.1.threads-$inx.txt
done
