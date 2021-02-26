#!/bin/tcsh

module purge PrgEnv-intel
module load PrgEnv-gnu openmpi

set verbose
setenv OMP_NUM_THREADS 1

setenv DATA_DIRECTORY ./data
setenv RESTART_DIRECTORY ./restart
setenv INPUT_DIRECTORY ./
setenv SOURCE_DIRECTORY ./
mkdir -p $DATA_DIRECTORY
mkdir -p $RESTART_DIRECTORY
cp input.f90 $DATA_DIRECTORY/

setenv MPI_TYPE_MAX 65536
setenv MPI_REQUEST_MAX 65536

srun -n 32 ./3dh 
echo 'Done'
