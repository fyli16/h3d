# Define shell
SHELL = csh

all:
	cd src; make
	cp src/h3d .

clean:
	-rm data/*
	-rm restart/*
	-rm mesh_vertices.dat
	-rm sim_id.txt
	cd src; make clean;
	# -rm 3dh

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







