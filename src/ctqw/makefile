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
	pylapack = -llapack
	
	#sed 's/^\tuse/!\tuse/g;s/dbesjn/bessel_jn/g' libctqw_intel.f90 > libctqw_gcc.f90
else
	FC = ifort
	fcompiler=intelem
	switch = -O2 -assume bscc
	libName = _intel
	lapack = -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
	pylapack = -lmkl_lapack95_lp64 -lmkl_rt -lpthread
endif

progname = ctqw$(libName)
objects = FFT.o libctqw$(libName).o fileOps.o main.o

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
	cp $(progname) $(BIN_PATH)
	rm $(progname)

# Fortran source code
$(progname): $(objects)
	$(FC) $(objects) -o $(progname) $(switch) $(lapack) $(fftw) $(matExp)
FFT.o: FFT.f90
	$(FC) -c $(switch) FFT.f90
libctqw$(libName).o: FFT.o libctqw$(libName).f90
	$(FC) -c $(switch) libctqw$(libName).f90
fileOps.o: fileOps.f90
	$(FC) -c $(switch) fileOps.f90
main.o: FFT.o libctqw$(libName).o fileOps.o main.f90
	$(FC) -c $(switch) main.f90

#Python wrappers
python: libpyctqw$(libName).so
	
libpyctqw$(libName).so: libctqw$(libName).f90
	f2py $(pylapack) $(fftw) $(matExp) -m libpyctqw$(libName) -c --fcompiler=$(fcompiler) libctqw$(libName).f90
	cp libpyctqw$(libName).so $(LIB_PATH)
	rm -f *.so

#Cleaning files
clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.so













