#!/usr/bin/make -f

#Makefile
all: libctqw pyctqw examples
fortran: libctqw examples
python: pyctqw

.PHONY: clean

include $(SLEPC_DIR)/conf/slepc_common
PETSC_INCLUDE = -I$(PETSC_DIR)/include
PETSC_ARCH_INCLUDE = -I$(PETSC_DIR)/$(PETSC_ARCH)/include

examples: bin/exampleMPI

bin/exampleMPI: lib/libctqwMPI.so examples/exampleMPI.o
	mkdir -p bin
	cd examples && $(FLINKER) exampleMPI.o -o ../bin/exampleMPI -L../lib -lctqwMPI $(SLEPC_LIB)

examples/exampleMPI.o: examples/exampleMPI.F90
	cd examples && $(FLINKER) -E exampleMPI.F90 $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) > temp.f90
	sed -i '/implicit none/,/! declare variables/d' examples/temp.f90
	sed -i '/# 1/,/!~/d' examples/temp.f90
	cd examples && $(FLINKER) -c temp.f90 -o exampleMPI.o -I../include
	rm examples/temp.f90
	
libctqw: lib/libctqwMPI.so
	
lib/libctqwMPI.so: src/libctqw-MPI.F90
#	cd src && $(FLINKER) -c libctqw-MPI.F90 -o libctqw-MPI.o $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE)
	#ar cr lib/libctqwMPI.a src/libctqw-MPI.o
	mkdir -p include
	mkdir -p lib
	cd src && $(FLINKER) -shared -fPIC -J../include libctqw-MPI.F90 -o ../lib/libctqwMPI.so $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE)

#Python wrappers	
pyctqw: src/libctqw-MPI.F90
	CC=${CC} F90=${FC} LDSHARED=${FC_LINKER} \
	python setup.py config_fc --f90flags="-Wno-unused-variable" -q install --user
	${RM} -rf src/libpyctqw_MPImodule.c src/libpyctqw_MPI-f2pywrappers2.f90 src/*.o src/*.mod
	
clean::
	rm -f */*.o */*.so */*.s */*.mod */*temp.f90
	rm -rf build lib include bin
	rm -rf src/libpyctqw_MPImodule.c src/libpyctqw_MPI-f2pywrappers2.f90

