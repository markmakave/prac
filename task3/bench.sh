#!/bin/bash

modes=("ijk" "ikj" "kij" "jik" "jki" "kji")

for m in ${modes[@]}
do
./main A.bin B.bin C.bin $m
done
