#!/bin/bash

for ratio in 100 120 140 160 180 200 250; do
    for value in 20 30 40 50 60 70 80 90 100; do
        for i in {1..1000}; do
            echo "Run $i with value $value,ratio $ratio"
            ./rand-dag $value $ratio && \
            timeout --foreground 10s /home/panagiotis/code/dag-partitioning/cmake-build-debug/dag_partitioning || \
            {
                ret=$?
                echo "Error: dag-partitioning failed for value $value, ratio $ratio on iteration $i"
                if [ $ret -eq 124 ]; then
                    echo "TIMED OUT"
                fi
                exit 1
            }
        done
    done
done