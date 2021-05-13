#!/bin/tcsh

# module purge
# set HOST=`uname -n`
# echo The host is $HOST.
# if ( $HOST =~ *wf-fe* ) then
#   module load gcc openmpi python-epd idl
#   else if ( $HOST =~ *gr-fe* ) then 
#     module load gcc/9.3.0 openmpi/3.1.6
#   else if ( $HOST =~ *frontera* ) then 
#     module load gcc/9.1.0 impi/19.0.9
#   else
#     echo "Unrecognize machine; nothing is done"
#   endif
# endif

module purge PrgEnv-intel
module load PrgEnv-gnu openmpi

echo "#--------------------- compling --------------------#"
make

# echo ""
# echo "#--------------------- cleaning --------------------#"
# make clean
