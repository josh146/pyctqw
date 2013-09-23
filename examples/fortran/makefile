#!/usr/bin/make -f

# If they are not automatically found, make sure
# your PETSc and SLEPc variables are correctly
# set by uncommenting out the below.
# are correctly
#PETSC_DIR = 
#PETSC_ARCH =
#SLEPC_DIR = 

CTQW_DIR = ../..
include $(CTQW_DIR)/ctqw_common

program = exampleMPI

#Makefile	
all: $(program)

$(program): $(program).o
	$(FLINKER) $(program).o -o $(program) $(CTQW_LIBS)

$(program).o: $(program).F90
	$(FLINKER) $(program).F90 $(CTQW_INCLUDE_DIRS) -c -o $(program).o 

clean::
	rm -f *.o *.mod *.s

.DEFAULT_GOAL := all