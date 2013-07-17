#!/usr/bin/make -f

progname = ctqw
objects = FFT.o libctqw.o fileOps.o main.o
BIN_INSTALL = bin
BIN_PATH = ../../$(BIN_INSTALL)
LIB_INSTALL = lib
LIB_PATH = ../../$(LIB_INSTALL)

ifeq ($(FC),gfortran)
	fcompiler=gfortran
	switch = -O2 -fbackslash
	libName = _gcc
	lapack = -llapack
	pylapack = -llapack
else
	FC = ifort
	fcompiler=intelem
	switch = -O2 -assume bscc
	libName = _intel
	lapack = -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
	pylapack = -lmkl_lapack95_lp64 -lmkl_rt -lpthread
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
libctqw.o: FFT.o libctqw.f90
	$(FC) -c $(switch) libctqw.f90
fileOps.o: fileOps.f90
	$(FC) -c $(switch) fileOps.f90
main.o: FFT.o libctqw.o fileOps.o main.f90
	$(FC) -c $(switch) main.f90

#Python wrappers
python: libpyctqw.so
	
libpyctqw.so: libctqw.f90
	f2py $(pylapack) $(fftw) $(matExp) -m libpyctqw -c --fcompiler=$(fcompiler) libctqw.f90
	cp libpyctqw.so $(LIB_PATH)

#Cleaning files
clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.so













