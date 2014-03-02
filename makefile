#!/usr/bin/make -f
#Makefile

CTQW_DIR = .

include $(CTQW_DIR)/ctqw_common

F90library = libctqwMPI

ifeq ($(shared_lib),1)
	EXT = so
else
	EXT = a
endif

ifeq ($(COMPILER),intel)
    MOD = -module
else
	MOD = -J
endif

# include
all: fortran examples

fortran: $(F90library)
examples: fexampleMPI

# Fortran MPI example
fexampleMPI:
	make -C $(CTQW_DIR)/examples/fortran
clean-fexampleMPI:
	make -C $(CTQW_DIR)/examples/fortran clean

# Fortran MPI library and interface
libctqwMPI: $(CTQW_LIB)/libctqwMPI.$(EXT)
	
lib/libctqwMPI.so: $(CTQW_DIR)/src/ctqwMPI.F90 
	mkdir -p $(CTQW_INCLUDE)
	mkdir -p $(CTQW_LIB)
	cd $(CTQW_DIR)/src && $(FLINKER) -shared -fPIC -Wno-conversion $(MOD) ../$(CTQW_INCLUDE) ctqwMPI.F90 -o ../$(CTQW_LIB)/libctqwMPI.so $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE)

lib/libctqwMPI.a: $(CTQW_DIR)/src/ctqwMPI.F90 
	mkdir -p $(CTQW_INCLUDE)
	mkdir -p $(CTQW_LIB)
	cd $(CTQW_DIR)/src && $(FLINKER) -Wno-conversion $(MOD) ../$(CTQW_INCLUDE) -c ctqwMPI.F90 -o libctqw-MPI.o $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE)
	ar cr $(CTQW_DIR)/lib/libctqwMPI.a $(CTQW_DIR)/src/libctqw-MPI.o

# documentation
.PHONY: docs
docs:
	make -C docs/ html
	ln -fs docs/_build/html/index.html ./index.html
docs-clean:
	make -C docs/ clean
	rm -f docs/stub/*

clean:: clean-fexampleMPI
	rm -f */*.o */*.so */*.s */*.mod
	rm -rf build lib include bin
	rm -rf src/libpyctqw_MPImodule.c src/libpyctqw_MPI-f2pywrappers2.f90
	rm -rf examples/fortran/exampleMPI
	rm -rf docs/_build/*
	rm -rf index.html

.DEFAULT_GOAL := all