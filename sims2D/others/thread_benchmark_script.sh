#!/bin/bash

# get the number of threads using the nproc command
NPROC=$(nproc)
for i in $(seq 1 $NPROC)
do
    echo "Running with $i threads"
    time julia --threads=$i thread_benchmark.jl
    echo ""
    rm -rf ./dump
done
