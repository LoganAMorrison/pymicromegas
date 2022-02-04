SHELL := /bin/bash

PYBIND11_INC = $(shell python3 -m pybind11 --includes)
EXTENSION = $(shell python3-config --extension-suffix)

.PHONY: all clean flags

all: include/microPath.h \
	build_calchep \
	build_micromegas \
	build_pymicromegas \
	move_pymicromegas

build_calchep: include/microPath.h
	$(MAKE) -C CalcHEP_src  MICROMEGAS=MICROMEGAS

build_micromegas: build_calchep
	$(MAKE) -C sources

build_pymicromegas: build_micromegas
	$(MAKE) -C MSSM --file=Makefile.pymicromegas

move_pymicromegas: build_pymicromegas
	# Move the libraries to the pymicromegas directory
	mv MSSM/pymicromegas$(EXTENSION) pymicromegas
	mv MSSM/micromegas$(EXTENSION) pymicromegas
	mv MSSM/softsusy$(EXTENSION) pymicromegas
	mv MSSM/spheno$(EXTENSION) pymicromegas
	mv MSSM/suspect$(EXTENSION) pymicromegas
	# Pre-moved the files
	# mv MSSM/pymicromegas/micromegas.pyi pymicromegas
	# mv MSSM/pymicromegas/pymicromegas.pyi pymicromegas
	# mv MSSM/pymicromegas/softsusy.pyi pymicromegas
	# mv MSSM/pymicromegas/suspect.pyi pymicromegas

install:
	$(shell pip install .)

include/microPath.h:
	echo \#define micrO \"$(CURDIR)\"  > include/microPath.h
      
   
clean:  
	./clean
flags: 
	$(MAKE) -C CalcHEP_src flags
