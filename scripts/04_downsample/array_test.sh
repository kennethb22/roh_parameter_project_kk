#!/bin/bash

declare -a cvgX=(50x 30x 15x 10x 05x)
declare -a cvgP=(1.0 0.6 0.3 0.2 0.1)
cvgCnt=${#cvgX[@]}
let cvgCnt-=1


for population in 100 50 30
do
    for i in $(seq 0 $cvgCnt)
    do 
        echo pop_${population}_cvg_${cvgP[i]}_cvgx_${cvgX[i]}
    done
done