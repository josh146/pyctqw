#!/usr/bin/make -f

progname = ctqw
objects = FFT.o myFunctions.o fileOps.o main.o
BIN_INSTALL = bin
BIN_PATH = ../../$(BIN_INSTALL)

ifeq ($(FC),gfortran)
	fcompiler=gfortran
	switch = -O2 -fbackslash
	libName = _gcc
	lapack = -llapack
else
	FC = ifort
	fcompiler=intelem
	switch = -O2 -assume bscc
	libName = _intel
	lapack = -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
endif

# fftw3 library location
ifeq ($(target),ivec)
	fftw = -L/ivec/devel/intel/12.1/fftw/3.3.3/lib -lfftw3
else
	fftw = -lfftw3
endif

# burkadt matrix exp library location
matExp = -L../../lib -lmatExp$(libName) -lc8$(libName) -lr8$(libName)

#Makefile
all: $(progname) clean
	mkdir -p $(BIN_PATH)
	cp ctqw $(BIN_PATH)
	rm ctqw

# Fortran source code
$(progname): $(objects)
	$(FC) $(objects) -o $(progname) $(switch) $(lapack) $(fftw) $(matExp)
FFT.o: FFT.f90
	$(FC) -c $(switch) FFT.f90
myFunctions.o: FFT.o myFunctions.f90
	$(FC) -c $(switch) myFunctions.f90
fileOps.o: fileOps.f90
	$(FC) -c $(switch) fileOps.f90
main.o: FFT.o myFunctions.o fileOps.o main.f90
	$(FC) -c $(switch) main.f90

#Python wrappers
python: pyctqw.so
	
pyctqw.so: myFunctions.f90
	f2py $(lapack) $(fftw) $(matExp) -m pyctqw -c --fcompiler=$(fcompiler) myFunctions.f90

#Cleaning files
clean:
	rm -f *.o
	rm -f *.mod













