# Define shell
SHELL = csh

all:
	cd src; make
	cp src/h3d .

clean:
	cd src; make clean

cleanall: clean
	cd src; make clean; make cleanall
	-rm data/hist*
	-rm data/*.gda
	-rm data/dat/*.gda
	-rm data/res*
	-rm restart_files/res*
	-rm output-*
	-rm 3dhout







