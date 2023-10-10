#!/bin/bash

modes=(0, 1, 2, 3, 4, 5)

for m in ${modes[@]}
do
./main A.bin B.bin C.bin $m
done
