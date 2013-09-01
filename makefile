#!/usr/bin/make -f

ifneq (${use_mpi},'False')
	include $(SLEPC_DIR)/conf/slepc_common
	PETSC_INCLUDE = -I$(PETSC_DIR)/include
	PETSC_ARCH_INCLUDE = -I$(PETSC_DIR)/$(PETSC_ARCH)/include
	use_mpi = True

	ex = F90examplesMPI
	F90library = libctqwMPI
else
	ex = F90examples
	F90library = libctqw
endif

#Makefile
all: $(F90library) pyctqw $(ex)

fortran: $(F90library)
examples: $(ex)
python: pyctqw

.PHONY: clean

# Fortran MPI example
F90examplesMPI: bin/exampleMPI

bin/exampleMPI: lib/libctqwMPI.so examples/exampleMPI.o
	mkdir -p bin
	cd examples && $(FLINKER) exampleMPI.o -o ../bin/exampleMPI -L../lib -lctqwMPI -lrefsor $(SLEPC_LIB)

examples/exampleMPI.o: examples/exampleMPI.F90
	cd examples && $(FLINKER) -E exampleMPI.F90 $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) > temp.f90
	sed -i '/implicit none/,/! declare variables/d' examples/temp.f90
	sed -i '/# 1/,/!~/d' examples/temp.f90
	cd examples && $(FLINKER) -c temp.f90 -o exampleMPI.o -I../include
	rm examples/temp.f90

# Fortran MPI library and interface
libctqwMPI: lib/libctqwMPI.so
	
lib/libctqwMPI.so: src/libctqw-MPI.F90 lib/librefsor.a
#	cd src && $(FLINKER) -c libctqw-MPI.F90 -o libctqw-MPI.o $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE)
	#ar cr lib/libctqwMPI.a src/libctqw-MPI.o
	mkdir -p include
	mkdir -p lib
	cd src && $(FLINKER) -shared -fPIC -Wno-conversion -J../include libctqw-MPI.F90 -o ../lib/libctqwMPI.so $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE) -I../include -L../include

lib/librefsor.a: src/sort.f90
	mkdir -p include
	mkdir -p lib
	cd src && $(FLINKER) -c sort.f90 -o sort.o
	ar cr lib/librefsor.a src/sort.o
	mv src/m_refsor.mod include/

#Python wrappers
pyctqw: src/libctqw-MPI.F90
	CC=${CC} F90=${FC} LDSHARED=${FC_LINKER}\
		use_mpi=$(use_mpi) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=$(PETSC_ARCH) SLEPC_DIR=$(SLEPC_DIR)\
		python setup.py config_fc --f90flags="-Wno-unused-variable" -q build

pyctqw\\ install: src/libctqw-MPI.F90
		python setup.py install $(prefix)
	${RM} -rf src/libpyctqw_MPImodule.c src/libpyctqw_MPI-f2pywrappers2.f90 src/*.o src/*.mod

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

	
clean::
	rm -f */*.o */*.so */*.s */*.mod */*temp.f90
	rm -rf build lib include bin
	rm -rf src/libpyctqw_MPImodule.c src/libpyctqw_MPI-f2pywrappers2.f90
	rm -rf docs/_build/*
	rm -rf index.html


