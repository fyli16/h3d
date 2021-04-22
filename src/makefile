# Choose compilers and options 
F90         = mpif90
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
OBJS =  parameter.of90 \
		utils.of90\
		c_int_to_char.ogcc \
		functions.of90 \
		mesh.of90 \
		field.of90\
		particle.of90 \
		eta_cal.of90 \
		boundary.of90 \
	 	diagnostics.of90\
		io.of90\
		init.of90 \
		h3d.of90

$(EFILE): $(OBJS)
	@echo "linking..."
	$(F90) $(OPTFLAGS2) -o $(EFILE) $(OBJS) $(LIBS)

# Compile C files 
%.ogcc: %.c 
	$(CC) -o $@ -c $<

# Compile F90 files
%.of90: %.f90
	$(F90) $(OPTFLAGS) -o $@ -c $< 

clean:
	-rm  *.mod *.of90 *.ogcc