#!/bin/bash
> /home/panagiotis/code/dag-partitioning/test/trace.txt
> /home/panagiotis/code/dagP/trace.txt

counter=1
for ratio in 100 150 200 250 300 500; do
    for value in 20 50 100 150 200 300 400 600 1000 2000; do
        partitions=()
        for ((i=2; i<=5; i++)); do
            partitions+=($i)
        done
        power=8
        while [ $power -le $value ]; do
            partitions+=($power)
            power=$(( power * 2 ))
        done
        for p in "${partitions[@]}"; do
            for i in {1..20}; do
                echo "Counter $counter, Run $i, Size $value, Ratio $ratio, Partitions $p"
                echo "Counter $counter, Run $i, Size $value, Ratio $ratio, Partitions $p" >> /home/panagiotis/code/dag-partitioning/test/trace.txt
                echo "Counter $counter, Run $i, Size $value, Ratio $ratio, Partitions $p" >> /home/panagiotis/code/dagP/trace.txt
                /home/panagiotis/code/dag-partitioning/test/rand-dag $value $ratio 0 && \
                timeout --foreground 5m /home/panagiotis/code/dag-partitioning/cmake-build-debug/dag_partitioning $p /home/panagiotis/code/dag-partitioning/test/dag.dot >> /home/panagiotis/code/dag-partitioning/test/trace.txt || \
                {
                    ret=$?
                    echo "Error: dag-partitioning failed for value $value, ratio $ratio, partitions $p on iteration $i"
                    if [ $ret -eq 124 ]; then
                        echo "TIMED OUT"
                    fi
                    exit 1
                }
                timeout --foreground 5m /home/panagiotis/code/dagP/cmake-build-debug/dagP /home/panagiotis/code/dag-partitioning/test/dag.dot $p >> /home/panagiotis/code/dagP/trace.txt || \
                {
                    ret=$?
                    echo "Error: dagP failed for value $value, ratio $ratio, partitions $p on iteration $i"
                    if [ $ret -eq 124 ]; then
                        echo "TIMED OUT"
                    fi
                }
            done
            ((counter ++))
        done
    done
done
