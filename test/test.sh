#!/bin/bash

> /home/panagiotis/code/dag-partitioning/test/trace.txt
counter=1
for ratio in 100 120 140 160 180 200 250 300 400 500 600; do
    for value in 20 30 40 50 60 70 80 90 100 120 140 160 180 200 250 300 350 400 500 600 700; do
        for i in {1..1000}; do
            echo "Counter $counter Run $i with value $value, ratio $ratio"
            echo "Counter $counter Run $i with value $value, ratio $ratio" >> /home/panagiotis/code/dag-partitioning/test/trace.txt
            ./rand-dag $value $ratio && \
            timeout --foreground 10s /home/panagiotis/code/dag-partitioning/cmake-build-debug/dag_partitioning >> /home/panagiotis/code/dag-partitioning/test/trace.txt || \
            {
                ret=$?
                echo "Error: dag-partitioning failed for value $value, ratio $ratio on iteration $i"
                if [ $ret -eq 124 ]; then
                    echo "TIMED OUT"
                fi
                exit 1
            }
        done
        ((counter ++))
    done
done