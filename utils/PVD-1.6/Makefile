############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Compilation de VTK_mod.f90 et VTR_mod.f90
## 
## Auteur          : Jalel Chergui (LIMSI-CNRS) <Jalel.Chergui@limsi.fr>
## Créé le         : Tue Aug 21 16:03:56 2007
## Dern. mod. par  : Jalel Chergui (LIMSI-CNRS) <Jalel.Chergui@limsi.fr>
## Dern. mod. le   : Wed May  6 09:37:58 2009
###############################################################################

# To be changed according to flavors
SHELL   = /bin/bash
FC	= ifort
LD	= ifort
FC	= gfortran
LD	= gfortran
FC	= mpif90
LD	= mpif90
FFLAGS  = -O3 #-fconvert=little-endian
# No modifications after this line.

SRCMOD	= VTK_mod.f90 VTR_mod.f90
OBJMOD	= VTK_mod.o VTR_mod.o

.SUFFIXES: .f90 .o

all: $(OBJMOD)

conv: convert.o $(OBJMOD)
	$(FC) $(FFLAGS) -o $@ $< $(OBJMOD)

end: endian.o
	$(FC) $(FFLAGS) -o $@ $<
%.o : %.f90
	$(FC) $(FFLAGS) -c $<

clean:;	rm -f $(OBJMOD) *.mod conv end
