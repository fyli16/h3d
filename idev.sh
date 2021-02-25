#!/bin/tcsh

module purge
module load gcc/9.1.0 impi/19.0.9

set verbose
setenv OMP_NUM_THREADS 1

setenv DATA_DIRECTORY ./data
setenv RESTART_DIRECTORY ./restart
setenv INPUT_DIRECTORY ./
setenv SOURCE_DIRECTORY ./
mkdir -p $DATA_DIRECTORY
mkdir -p $RESTART_DIRECTORY
cp input.f90 $DATA_DIRECTORY/

ibrun -np 56 $SOURCE_DIRECTORY/3dh 
echo 'Done'
