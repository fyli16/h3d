#!/bin/tcsh

module purge
module load impi

# set verbose
setenv OMP_NUM_THREADS 1

mkdir -p data; mkdir -p restart

mkdir -p data/bx;  mkdir -p data/by;  mkdir -p data/bz
mkdir -p data/ex;  mkdir -p data/ey;  mkdir -p data/ez
mkdir -p data/den; mkdir -p data/eta; mkdir -p data/eta_par

mkdir -p data/particle; mkdir -p data/probes

mkdir -p data/fox;  mkdir -p data/foy;   mkdir -p data/foz
mkdir -p data/p-xx; mkdir -p data/p-xy;  mkdir -p data/p-xz
mkdir -p data/p-yy; mkdir -p data/p-yz;  mkdir -p data/p-zz
mkdir -p data/tpar; mkdir -p data/tperp; mkdir -p data/tracking
mkdir -p data/vix;  mkdir -p data/viy;   mkdir -p data/viz 
mkdir -p data/vxs;  mkdir -p data/vys;   mkdir -p data/vzs 

mkdir -p data/ecal

# setenv MPI_TYPE_MAX 65536
# setenv MPI_REQUEST_MAX 65536

srun -n 32 ./build/h3d
exit