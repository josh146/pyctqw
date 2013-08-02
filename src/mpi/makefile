#!/usr/bin/make -f

#Makefile
all: example python

.PHONY: clean

include $(SLEPC_DIR)/conf/slepc_common
PETSC_INCLUDE = -I/opt/petsc/petsc-3.4.0/include

example: libctqw-MPI.o example.o
	$(FLINKER) libctqw-MPI.o example.o -o example $(SLEPC_LIB)

example.o: example.F90
	$(FLINKER) -E example.F90 $(PETSC_INCLUDE) > temp.f90
	sed -i '/implicit none/,/! declare variables/d' temp.f90
	sed -i '/# 1/,/!~/d' temp.f90
	$(FLINKER) -c temp.f90 -o example.o
	rm temp.f90
	
libctqw-MPI.o: libctqw-MPI.F90
	$(FLINKER) -c libctqw-MPI.F90 -o libctqw-MPI.o $(PETSC_INCLUDE) $(SLEPC_INCLUDE)

#Python wrappers	
python: libctqw-MPI.so

libctqw-MPI.so: libctqw-MPI.F90
	CC=${CC} F90=${FC} LDSHARED=${FC_LINKER} \
	python setup.py config_fc --f90flags="-Wno-unused-variable" -q build_ext --inplace
	${RM} -r build libpyctqw_MPImodule.c libpyctqw_MPI-f2pywrappers2.f90
	
clean::
	rm -f *.mod
	rm -f *.o *.so *.s
	rm -f *temp.f90 example
	rm -f *.py[co]
	rm -f -r __pycache__
	rm -f *.c *.f90
	rm -rf out/

