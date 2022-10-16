#!/bin/bash

# get the number of threads using the nproc command
NPROC=$(nproc)
for i in $(seq 1 $NPROC)
do
    echo "Running with $i threads" >> thread_benchmark.txt
    time julia.exe --threads=$i thread_benchmark.jl >> thread_benchmark.txt
    echo "" >> thread_benchmark.txt
    rm -rf ./dump
done