#!/bin/tcsh
#SBATCH -J h3d
#SBATCH -o log.%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 56
#SBATCH -t 01:00:00
#SBATCH -A PHY20020     

date

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

ibrun -np $SLURM_NTASKS ./h3d 
exit
