#!/bin/bash

MONITOR_PID=""
PID=""

cleanup() {
    echo -e "\n\nCleaning up..."
    
    if [[ -n "$MONITOR_PID" ]]; then
        kill "$MONITOR_PID" 2>/dev/null
        wait "$MONITOR_PID" 2>/dev/null
    fi
    
    if [[ -n "$PID" ]]; then
        kill "$PID" 2>/dev/null
        wait "$PID" 2>/dev/null
    fi

    
    printf "\033[?25h"
    rm -f trace.txt
    rm -f dag*dot
    rm -f dag*dot.node-mappings.txt
    
    echo "Cleanup completed. Exiting..."
    exit 130
}

trap cleanup SIGINT SIGTERM

clear
printf "\033[2J\033[H"

run_tests() {
    local minpar_list="$1"
    local output_file="$2"
    
    > "$output_file"
    
    local counter=0
    for minpar in $minpar_list; do
        for value in 20 50 100 200 500 1000 2000 5000 10000 50000 100000; do
            partitions=()
            power=2
            while [ $power -le $value ]; do
                partitions+=($power)
                power=$(( power * 2 ))
            done
            for p in "${partitions[@]}"; do
                for ratio in 100 120 150 200; do
                    for i in {1..5}; do
                        if [[ ! -f "$output_file" ]]; then
                            return 1
                        fi
                        
                        echo "Counter $counter, Run $i, Minpar $minpar, Size $value, Partitions $p, Ratio $ratio" >> "$output_file"
                        ../install/bin/rand-dag $value $ratio 0 ./dag-$counter_offset.dot && \
                        timeout --foreground 5m ../install/bin/test $p dag-$counter_offset.dot 1 $minpar >> "$output_file" || \
                        {
                            ret=$?
                            echo "Error: dag-partitioning failed for minpar $minpar, size $value, partitions $p, $ratio, on iteration $i"
                            if [ $ret -eq 124 ]; then
                                echo "TIMED OUT"
                            fi
                            exit 1
                        }
                    done
                    ((counter++))
                done
            done
        done
    done
}

monitor_progress() {
    local file="$1"
    local total="$2"
    
    printf "\033[?25l"
    
    while true; do
        printf "\033[H"
        
        if [[ -f $file ]]; then
            last_line=$(grep "^Counter" $file 2>/dev/null | tail -n 1)
            if [[ -n "$last_line" ]]; then
                current_counter=$(echo "$last_line" | sed -n 's/^Counter \([0-9]*\),.*/\1/p')

                progress=$((current_counter * 100 / total))
                bar_length=$((current_counter * 40 / total))
            else
                progress=0
                bar_length=0
            fi
        else
            progress=0
            bar_length=0
        fi
        printf "\033[KProgress: ["
        printf "%${bar_length}s" | tr ' ' '='
        printf "%$((40 - bar_length))s" | tr ' ' '.'
        printf "] %3d%% (%d/%d)\n" "$progress" "$current_counter" "$total"
        
        sleep 1
    done
}

> trace.txt

monitor_progress trace.txt 4200 &
MONITOR_PID=$!

(run_tests "0 10 20 50 100 200 500 1000 2000 5000" trace.txt 1) &
PID=$!

wait $PID
STATUS=$?

kill $MONITOR_PID 2>/dev/null
wait $MONITOR_PID 2>/dev/null
printf "\033[?25h\033[4B\n"

if [ $STATUS -ne 0 ]; then
    echo "Job failed with status $STATUS"
    rm -f trace.txt
    rm -f dag*dot
    rm -f dag*dot.node-mappings.txt
    exit 1
fi

rm -f dag*dot
rm -f dag*dot.node-mappings.txt

echo "All tests completed successfully"
