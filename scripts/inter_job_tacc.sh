#!/bin/tcsh
module purge
module load gcc/9.1.0 impi/19.0.9
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

ibrun -np 56 ./h3d
exit
