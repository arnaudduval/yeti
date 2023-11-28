#!/bin/bash

rm time.txt

. ~/py310yeti/bin/activate

for i in 1 2 4 6 7 8 12
do
    export OMP_NUM_THREADS=${i}
    python bench_solver.py
done