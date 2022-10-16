#!/bin/bash

# get the number of threads using the nproc command
NPROC=$(nproc)
for i in $(seq 1 $NPROC)
do
    echo "Running with $i threads" >> thread_benchmark.time
    { time julia.exe --threads=$i thread_benchmark.jl 2> thread_benchmark.out } 2> thread_benchmark.time
    echo "" >> thread_benchmark.time
    # rm thread_benchmark.out
    rm -rf ./dump
done
