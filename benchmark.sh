#!/bin/bash

PROGRAM="./integration"
FUNCTIONS=(1 2)
MODES=(1 2)
TOLERANCES=("1e-6" "1e-8")
PROCESSES=(1 2 4 8)
ITERATIONS=10

OUTPUT_FILE="benchmark_results.csv"
echo "Function,Mode,Tolerance,Processes,AverageTime(s)" > $OUTPUT_FILE

make clean && make

for func in "${FUNCTIONS[@]}"; do
    for mode in "${MODES[@]}"; do
        for tol in "${TOLERANCES[@]}"; do
            for p in "${PROCESSES[@]}"; do
                if [ "$mode" -eq 2 ] && [ "$p" -eq 1 ]; then
                    echo "Skipping Mode 2 for P=1 (requires at least 2 processes)"
                    continue
                fi
    

                echo "Testing: Func $func, Mode $mode, Tol $tol, P $p"
                total_time=0
                
                for ((i=1; i<=ITERATIONS; i++)); do                    
                    
                    raw_output=$(mpirun --oversubscribe -np $p $PROGRAM $func $mode $tol 2>&1)                                                        
                    exec_time=$(echo "$raw_output" | grep "Time:" | grep -oE '[0-9]+\.[0-9]+')
                    
                    if [ "$i" -eq 1 ]; then
                        echo "  Run 1 (Warm-up): ${exec_time:-FAILED} s"
                    else
                        if [[ -n "$exec_time" ]]; then
                            total_time=$(echo "$total_time + $exec_time" | bc -l)
                        else
                            echo "  Warning: Run $i failed to return a time."
                        fi
                    fi
                done
                
                average_time=$(echo "scale=6; $total_time / 9" | bc -l)
                echo "$func,$mode,$tol,$p,$average_time" >> $OUTPUT_FILE
                printf "  Average Time: %.6f s\n" $average_time
            done
        done
    done
done