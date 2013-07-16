#Makefile for CtQW
#Variables
progname = CtQW
mf = main
compiler = ifort
objects = FFT.o myFunctions.o fileOps.o $(mf).o
switch = -O2 -assume bscc

fftw = -L/ivec/devel/intel/12.1/fftw/3.3.3/lib -lfftw3
lapack = -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
matExp = -L/home/jizaac/lib -lmatExp_intel -lc8_intel -lr8_intel

#Makefile
$(progname): $(objects)
	$(compiler) $(objects) -o $(progname) $(switch) $(lapack) $(fftw) $(matExp)
FFT.o: FFT.f90
	$(compiler) -c $(switch) FFT.f90
myFunctions.o: FFT.o myFunctions.f90
	$(compiler) -c $(switch) myFunctions.f90
fileOps.o: fileOps.f90
	$(compiler) -c $(switch) fileOps.f90
$(mf).o: FFT.o myFunctions.o fileOps.o $(mf).f90
	$(compiler) -c $(switch) $(mf).f90
#Cleaning files
clean:
	rm *.o
	rm *.mod
	rm $(progname)
