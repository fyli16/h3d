##### Specify where the code will be run #####


######## Choose compilers and options ########

F90         = mpif90
CC          = mpicc
#OPTFLAGS    = -O3 
#OPTFLAGS2   = -O3
   OPTFLAGS    = -O3 -ffree-line-length-none -fimplicit-none -cpp
   OPTFLAGS2   = -O3 -ffree-line-length-none -fimplicit-none  
#   OPTFLAGS    = -Og -fbounds-check -ffree-line-length-none -fimplicit-none -g -cpp
#   OPTFLAGS2   = -Og -fbounds-check -ffree-line-length-none -fimplicit-none -g 
   LIBS        = 

############## Define dependecies ############

all: tracking picking

tracking: tracking.f90
	$(F90) $(OPTFLAGS2) -o $@ $< $(LIBS)

picking: picking.f90
	$(F90) $(OPTFLAGS2) -o $@ $< $(LIBS)


############## Compile C files ###############

%.ogcc: %.c makefile 
	$(CC) -o $@ -c $<

%.c: %.c.nashi makefile 
	cp -f $< $@

############## Compile F90 files #############

%.of90: %.f90 makefile 
	$(F90)  $(OPTFLAGS) -o $@ -c $< 

clean:
	-rm  *.mod  3dh *.of90 *.ogcc
cleanall: clean
	-rm  restart_data.*
   
