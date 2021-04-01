# Specify where the code will be run


# Choose compilers and options 
F90         = mpif90
CC          = mpicc
# OPTFLAGS    = -O3 
# OPTFLAGS2   = -O3
OPTFLAGS    = -O3 -ffree-line-length-none -fimplicit-none -cpp
OPTFLAGS2   = -O3 -ffree-line-length-none -fimplicit-none  
# OPTFLAGS    = -Og -fbounds-check -ffree-line-length-none -fimplicit-none -g -cpp
# OPTFLAGS2   = -Og -fbounds-check -ffree-line-length-none -fimplicit-none -g 
LIBS        = 

C_FIELS = ./C

# Define dependecies
EFILE = h3d
OBJS =  parameter.of90 \
		mesh_class.of90 \
		mesh2d.of90 \
		functions.of90 \
		init_single_alfven_wave.of90 \
		init_ion_acoustic_wave.of90 \
		init_lapd_alfven_wave.of90 \
		io.of90\
		utils.of90\
	 	boundary.of90 \
		field.of90\
		wrtfile.of90 \
		wrtfile_Non_MPIO.of90\
		user_diagnostics.of90\
		initialize.of90\
		3dHybrid.of90 \
		${C_FIELS}/get_env_var.ogcc \
		${C_FIELS}/integer_to_character.ogcc \
		particle.of90 \
		field_2d.of90 \
		etacal.of90 \
		etacal_2d.of90\
		inject_wave.of90

$(EFILE): $(OBJS)
	@echo "linking..."
	$(F90) $(OPTFLAGS2) -o $(EFILE) $(OBJS) $(LIBS)

# Compile C files 
%.ogcc: %.c 
	$(CC) -o ${C_FIELS}/$@ -c ${C_FIELS}/$<

# Compile F90 files
%.of90: %.f90
	$(F90) $(OPTFLAGS) -o $@ -c $< 

clean:
	-rm  *.mod *.of90 ${C_FIELS}/*.ogcc h3d
cleanall: clean
	-rm  restart_data.*