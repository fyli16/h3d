#!/bin/tcsh

module purge PrgEnv-intel
module load PrgEnv-gnu openmpi
echo ""
echo "#--------------------- compling --------------------#"
make
echo ""
echo "#--------------------- cleaning --------------------#"
make clean
echo ""


set verbose
setenv OMP_NUM_THREADS 1

setenv DATA_DIRECTORY ./data
setenv RESTART_DIRECTORY ./restart
mkdir -p $DATA_DIRECTORY
mkdir -p $RESTART_DIRECTORY

# setenv MPI_TYPE_MAX 65536
# setenv MPI_REQUEST_MAX 65536

srun -n 32 ./src/h3d
exit