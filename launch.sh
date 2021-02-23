#!/bin/sh

for i in 0.001 0.01 0.05 0.1 0.2 0.5
do 
    echo "for b=${i} ..."
    mkdir ../1d-b${i}
    cd ../1d-b${i}
    cp -r ../launch/src .
    cp ../launch/Makefile .
    cp ../launch/finput.f90 .
    cp ../launch/job_tacc .
    sed -i "s/dB_B0=0.1/dB_B0=${i}/" src/init_single_alfven_wave.f90
    sbatch job_tacc
    cd ../launch
done