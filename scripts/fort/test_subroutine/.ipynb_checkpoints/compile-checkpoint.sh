gfortran -o param.of90 -c param.f90
gfortran -o sub1.of90 -c sub1.f90
gfortran test_subroutine.f90 sub1.of90 param.of90

# rm *.of90
