#!/bin/bash

for ratio in 100 150 200 300 500 750 1000; do
    for value in 20 50 100 200 500 1000 2000 5000; do
        for i in {1..200}; do
            echo "Run $i with value $value,ratio $ratio"
            ./rand-dag $value $ratio && \
            /home/panagiotis/code/dag-partitioning/cmake-build-debug/dag_partitioning || \
            {
                echo "Error: dag-partitioning failed for value $value, ratio $ratio on iteration $i"
                exit 1
            }
        done
    done
done