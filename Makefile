# Define shell
SHELL = csh

all:
	cd src; make
	cp src/3dh .

clean:
	cd src; make clean;
	# -rm 3dh
	# -rm data/*
	# -rm restart_files/*
	# -rm mesh_vertices.dat

cleanall: clean
	cd src; make clean; make cleanall
	# cd post-process; make clean 
	-rm data/hist*
	-rm data/*.gda
	-rm data/dat/*.gda
	-rm data/res*
	-rm restart_files/res*
	-rm output-*
	-rm 3dhout







