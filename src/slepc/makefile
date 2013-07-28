#!/usr/bin/make -f

#Makefile
all: slepc

.PHONY: clean

include ${SLEPC_DIR}/conf/slepc_common
PETSC_INCLUDE = -I/opt/petsc/petsc-3.4.0/include

slepc: slepc.o
	${FLINKER} slepc.o -o slepc ${SLEPC_LIB}

slepc.o: slepc.F90
	${FLINKER} -c slepc.F90 -o slepc.o ${PETSC_INCLUDE} ${SLEPC_INCLUDE} ${SLEPC_LIB}
	
clean::
	rm -f *.mod
	rm -f slepc

