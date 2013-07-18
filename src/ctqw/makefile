#!/usr/bin/make -f

BIN_INSTALL = bin
BIN_PATH = ../../$(BIN_INSTALL)
LIB_INSTALL = lib
LIB_PATH = ../../$(LIB_INSTALL)

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
	lapack = -lmkl_lapack95_lp64 -lmkl_rt -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
endif

progname = $(BIN_PATH)/ctqw$(libName)
pylibpath = $(LIB_PATH)/libpyctqw$(libName).so
objects = FFT$(libName).o libctqw$(libName).o fileOps$(libName).o main$(libName).o

# fftw3 library location
ifeq ($(target),ivec)
	fftw = -L/ivec/devel/intel/12.1/fftw/3.3.3/lib -lfftw3
else
	fftw = -lfftw3
endif

# burkadt matrix exp library location
matExp = -L../../lib -lmatExp$(libName) -lc8$(libName) -lr8$(libName)

#Makefile
all: $(progname)

# Fortran source code
$(progname): $(objects)
	$(FC) $(objects) -o $(progname) $(switch) $(lapack) $(fftw) $(matExp)
	rm -f *.mod
FFT$(libName).o: FFT.f90
	$(FC) -c $(switch) FFT.f90 -o FFT$(libName).o
libctqw$(libName).o: FFT$(libName).o libctqw$(libName).f90
	$(FC) -c $(switch) libctqw$(libName).f90 -o libctqw$(libName).o
fileOps$(libName).o: fileOps.f90
	$(FC) -c $(switch) fileOps.f90 -o fileOps$(libName).o
main$(libName).o: FFT$(libName).o libctqw$(libName).o fileOps$(libName).o main.f90
	$(FC) -c $(switch) main.f90 -o main$(libName).o
	
libctqw_gcc.f90: libctqw_intel.f90
	echo "making gcc fortran file"
	sed 's/^\tuse/!\tuse/g;s/dbesjn/bessel_jn/g;/^\tsubroutine progressbar.*/,/^\tend subroutine progressbar.*/s/^/!/;/call progressbar(m,terms)/s/^/!/' libctqw_intel.f90 > libctqw_gcc.f90

#Python wrappers
python: $(pylibpath)
	
$(pylibpath): libctqw$(libName).f90
	f2py $(lapack) $(fftw) $(matExp) -m libpyctqw$(libName) -c --fcompiler=$(fcompiler) libctqw$(libName).f90
	cp libpyctqw$(libName).so $(LIB_PATH)
	rm -f *.so
	rm -f *.mod

#Cleaning files
clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.so













