#!/usr/bin/make -f

#Makefile
all: ctqwMPI python

.PHONY: clean

include $(SLEPC_DIR)/conf/slepc_common
PETSC_INCLUDE = -I/opt/petsc/petsc-3.4.0/include

ctqwMPI: libctqwMPI.o ctqwMPI.o
	$(FLINKER) libctqwMPI.o ctqwMPI.o -o ctqwMPI $(SLEPC_LIB)

ctqwMPI.o: ctqwMPI.F90
	$(FLINKER) -E ctqwMPI.F90 $(PETSC_INCLUDE) > temp.f90
	sed -i '/implicit none/,/! declare variables/d' temp.f90
	sed -i '/# 1/,/!~/d' temp.f90
	$(FLINKER) -c temp.f90 -o ctqwMPI.o
	rm temp.f90
	
libctqwMPI.o: libctqwMPI.F90
	$(FLINKER) -c libctqwMPI.F90 -o libctqwMPI.o $(PETSC_INCLUDE) $(SLEPC_INCLUDE)

#Python wrappers
FILTER_OUT = $(foreach v,$(2),$(if $(findstring $(1),$(v)),,$(v)))

F90FLAGS = $(call FILTER_OUT,mpif90, $(FLINKER))
F90FLAGS += $(SLEPC_LIB) $(PETSC_INCLUDE) $(SLEPC_INCLUDE)

F2PY_FLAGS = --quiet --noarch -DF2PY_REPORT_ON_ARRAY_COPY=1
F2PY_FLAGS += --fcompiler=gnu95 --f90exec=mpif90 --f90flags="$(F90FLAGS)"
	
python: libctqwMPI.so

libctqwMPI.so: libctqwMPI.F90
	CC=${CC} F90=${FC} LDSHARED=${FC_LINKER} \
	python setup.py config_fc --f90flags="-Wno-unused-variable" -q build_ext --inplace
	${RM} -r build libctqwMPImodule.c libctqwMPI-f2pywrappers2.f90
	
clean::
	rm -f *.mod
	rm -f *.o *.so *.s
	rm -f ctqwMPI *temp.f90
	rm -f *.py[co]
	rm -f -r __pycache__
	rm -f *.c *.f90
	rm -rf out/

