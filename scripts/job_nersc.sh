#!/bin/tcsh

#SBATCH -J h3d      
#SBATCH -o log.%j
#SBATCH --qos regular       
#SBATCH -N 2            
#SBATCH -n 64       
#SBATCH -t 01:30:00 
#SBATCH -A m2407  #MR
##SBATCH -A m3757  #CSR
#SBATCH --constraint=haswell

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

srun ./src/h3d
exit