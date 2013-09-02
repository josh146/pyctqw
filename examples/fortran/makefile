#!/usr/bin/make -f

# If they are not automatically found, make sure
# your PETSc and SLEPc variables are correctly
# set by uncommenting out the below.
# are correctly
#PETSC_DIR = 
#PETSC_ARCH =
#SLEPC_DIR = 

CTQW_DIR = ../..

CTQW_LIB = $(CTQW_DIR)/lib
CTQW_INCLUDE = $(CTQW_DIR)/include

program = exampleMPI

#Makefile	
all: $(program)

# include the CTQW directives
include $(CTQW_DIR)/ctqw_common

$(program): $(program).o
	$(FLINKER) $(program).o -o $(program) -L$(CTQW_LIB) -lctqwMPI $(SLEPC_LIB)

# Preprocess the F90 Fortran program, removes duplicate declarations, and compiles
# IMPORTANT: node the placement of CTQW_DIRECTIVES, between the source file and '-c -o'.
# 	THIS ORDERING IS IMPORTANT.
$(program).o: $(program).F90
	$(FLINKER) $(program).F90 $(CTQW_DIRECTIVES) -c -o $(program).o -I$(CTQW_INCLUDE)

clean::
	rm -f *.o *.mod *.s

.DEFAULT_GOAL := all