F90 = mpif90
EFILE = test_mpi
OPTFLAGS = -O3 -ffree-line-length-none -fimplicit-none 
OBJS = cal_pi.of90

$(EFILE): $(OBJS)
	@echo "linking ..."
	$(F90) $(OPTFLAGS) -o $(EFILE) $(OBJS)

%.of90: %.f90 
	$(F90) $(OPTFLAGS) -o $@ -c $< 

clean:
	-rm *.of90