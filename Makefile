SHELL := /bin/bash

PYBIND11_INC = $(shell python3 -m pybind11 --includes)
EXTENSION = $(shell python3-config --extension-suffix)

.PHONY: all

all: build_micromegas \
	build_pymicromegas \
	move_pymicromegas

build_micromegas:
	$(MAKE) -C micromegas

build_pymicromegas: build_micromegas
	$(MAKE) -C micromegas/MSSM --file=Makefile.pymicromegas

move_pymicromegas: build_pymicromegas
	# Move the libraries to the pymicromegas directory
	mv micromegas/MSSM/pymicromegas$(EXTENSION) pymicromegas
	mv micromegas/MSSM/micromegas$(EXTENSION) pymicromegas
	mv micromegas/MSSM/softsusy$(EXTENSION) pymicromegas
	mv micromegas/MSSM/spheno$(EXTENSION) pymicromegas
	mv micromegas/MSSM/suspect$(EXTENSION) pymicromegas

install:
	$(shell pip install .)
