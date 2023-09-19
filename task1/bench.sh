#!/bin/bash

threads=(1, 2, 3, 4)

for t in ${threads[@]}
do
./main 100000000 $t
done
