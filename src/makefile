# Choose compilers and options 
FC         = mpif90
CC          = mpicc
# OPTFLAGS    = -O3 
# OPTFLAGS2   = -O3
OPTFLAGS    = -O3 -ffree-line-length-none -fimplicit-none -cpp #-fcheck=all 
OPTFLAGS2   = -O3 -ffree-line-length-none -fimplicit-none #-fcheck=all 
# OPTFLAGS    = -Og -fbounds-check -ffree-line-length-none -fimplicit-none -g -cpp
# OPTFLAGS2   = -Og -fbounds-check -ffree-line-length-none -fimplicit-none -g 
LIBS        = 

# Define dependecies
EFILE = h3d
OBJS =  param.of90 \
				utils.of90\
				int2char.ogcc \
				func.of90 \
				mesh.of90 \
				field.of90\
				particle.of90 \
				eta.of90 \
				boundary.of90 \
				io.of90\
				restart.of90 \
				diag.of90\
				init.of90 \
				h3d.of90

$(EFILE): $(OBJS)
	@echo "linking..."
	$(FC) $(OPTFLAGS2) -o $(EFILE) $(OBJS) $(LIBS)

# Compile C files 
%.ogcc: %.c 
	$(CC) -o $@ -c $<

# Compile Fortran files
%.of90: %.f90
	$(FC) $(OPTFLAGS) -o $@ -c $< 

clean:
	-rm  *.mod *.of90 *.o *.ogcc