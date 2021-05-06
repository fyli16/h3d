#!/bin/sh
module purge PrgEnv-intel
module load PrgEnv-gnu openmpi

gfortran -o constants.of90 -c constants.f90
gfortran -o cal_area.of90 -c cal_area.f90
gfortran main.f90 constants.of90 cal_area.of90

# rm *.of90 *.mod
