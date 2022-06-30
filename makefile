#!/bin/tcsh

# Choose compilers and options 
FC          = mpif90
CC          = mpicc

OPTFLAGS    = -O3 -ffree-line-length-none -fimplicit-none -cpp #-fcheck=all 
OPTFLAGS2   = -O3 -ffree-line-length-none -fimplicit-none #-fcheck=all 

LIBS        = 

BUILD_DIR = ./build
$(shell mkdir -p $(BUILD_DIR))
SRC_DIR   = ./src

# Define dependecies
EFILE = h3d
OBJS =  $(BUILD_DIR)/param.of90 \
				$(BUILD_DIR)/utils.of90\
				$(BUILD_DIR)/int2char.ogcc \
				$(BUILD_DIR)/func.of90 \
				$(BUILD_DIR)/mesh.of90 \
				$(BUILD_DIR)/injection.of90\
				$(BUILD_DIR)/field.of90\
				$(BUILD_DIR)/particle.of90 \
				$(BUILD_DIR)/eta.of90 \
				$(BUILD_DIR)/boundary.of90 \
				$(BUILD_DIR)/io.of90\
				$(BUILD_DIR)/restart.of90 \
				$(BUILD_DIR)/diag.of90\
				$(BUILD_DIR)/init.of90 \
				$(BUILD_DIR)/h3d.of90

$(EFILE): $(OBJS)
	@echo "linking..."
	$(FC) $(OPTFLAGS2) -o $(BUILD_DIR)/$(EFILE) $(OBJS) $(LIBS)

# Compile C files 
$(BUILD_DIR)/%.ogcc: $(SRC_DIR)/%.c 
	$(CC) -o $@ -c $<

# Compile Fortran files
$(BUILD_DIR)/%.of90: $(SRC_DIR)/%.f90
	$(FC) $(OPTFLAGS) -J$(BUILD_DIR) -o $@ -c $< 

clean:
	-rm $(BUILD_DIR)/*.mod 
	-rm $(BUILD_DIR)/*.of90
	-rm $(BUILD_DIR)/*.o 
	-rm $(BUILD_DIR)/*.ogcc 
	-rm $(BUILD_DIR)/${EFILE}







