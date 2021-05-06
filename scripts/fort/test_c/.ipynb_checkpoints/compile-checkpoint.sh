#!/bin/sh
module purge PrgEnv-intel
module load PrgEnv-gnu openmpi

# mpicc -o int_to_char.ogcc -c int_to_char.c
# mpif90 -O3 -ffree-line-length-none -fimplicit-none -cpp -o test_c.of90 -c test_c.f90
# mpif90 -O3 -ffree-line-length-none -fimplicit-none -o test test_c.of90 int_to_char.ogcc

gcc -o int_to_char.ogcc -c int_to_char.c
gfortran -O3 -ffree-line-length-none -fimplicit-none -cpp -o test_c.of90 -c test_c.f90
gfortran -O3 -ffree-line-length-none -fimplicit-none -o test test_c.of90 int_to_char.ogcc

# rm *.of90 *.mod
