#!/usr/bin/make -f

# ifneq (${use_mpi},'False')
# 	include $(SLEPC_DIR)/conf/slepc_common
# 	PETSC_INCLUDE = -I$(PETSC_DIR)/include
# 	PETSC_ARCH_INCLUDE = -I$(PETSC_DIR)/$(PETSC_ARCH)/include
# 	use_mpi = True

# 	ex = fexampleMPI
# 	F90library = libctqwMPI
# else
# 	ex = fexample
# 	F90library = libctqw
# endif

#Makefile

CTQW_DIR = .

CTQW_LIB = $(CTQW_DIR)/lib
CTQW_INCLUDE = $(CTQW_DIR)/include

F90library = libctqwMPI

# include
all: fortran examples
 
include $(CTQW_DIR)/ctqw_common

fortran: $(F90library)
examples: fexampleMPI

# Fortran MPI example
fexampleMPI:
	make -C $(CTQW_DIR)/examples/fortran
clean-fexampleMPI:
	make -C $(CTQW_DIR)/examples/fortran clean

# Fortran MPI library and interface
libctqwMPI: $(CTQW_LIB)/libctqwMPI.so
	
lib/libctqwMPI.so: $(CTQW_DIR)/src/libctqw-MPI.F90 
	mkdir -p $(CTQW_INCLUDE)
	mkdir -p $(CTQW_LIB)
	cd $(CTQW_DIR)/src && $(FLINKER) -shared -fPIC -Wno-conversion -J../$(CTQW_INCLUDE) libctqw-MPI.F90 -o ../$(CTQW_LIB)/libctqwMPI.so $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE)

# documentation
docs-html:
	make -C docs/ html
	ln -fs docs/_build/html/index.html ./index.html
docs-pdf:
	make -C docs/ latex
	make -C docs/_build/latex
	mv docs/_build/latex/pyCTQW.pdf ./
clean-docs: clean
	make -C docs/ clean
	rm -f docs/stub/*

clean:: clean-fexampleMPI
	rm -f */*.o */*.so */*.s */*.mod
	rm -rf build lib include bin
	rm -rf src/libpyctqw_MPImodule.c src/libpyctqw_MPI-f2pywrappers2.f90
	rm -rf docs/_build/*
	rm -rf index.html
