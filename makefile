#!/usr/bin/make -f

BIN_INSTALL = bin
LIB_INSTALL = lib

ifeq ($(FC),gfortran)
	libName = _gcc
else
	FC = ifort
	libName = _intel
endif

progname = ctqw$(libName)
pylibpath = $(LIB_INSTALL)/libpyctqw$(libName).so
binary = $(BIN_INSTALL)/$(progname)

#sourcecode = $(shell ls src/ctqw/*.f90)
sourcecode = src/ctqw/FFT.f90 src/ctqw/libctqw$(libName).f90 src/ctqw/fileOps.f90 src/ctqw/main.f90
LibmatExp = $(LIB_INSTALL)/libr8$(libName).a $(LIB_INSTALL)/libc8$(libName).a $(LIB_INSTALL)/libmatExp$(libName).a

#Makefile
all: $(binary) $(pylibpath)

$(binary): $(LibmatExp) $(sourcecode)
	$(MAKE) -C src/ctqw
	
$(LibmatExp):
	$(MAKE) -C src/burkadt

$(sourcecode):
	
# Python library wrappers
libpyctqw: $(pylibpath)

$(pylibpath): src/ctqw/libctqw$(libName).f90 $(LibmatExp) 
	$(MAKE) python -C src/ctqw
	
# python CLI
pyctqw: $(pylibpath)
	$(MAKE) -C src/python

#Cleaning files
clean:
	$(MAKE) clean -C src/ctqw
	$(MAKE) clean -C src/python
	rm -f $(binary)

fullclean: clean
	rm -rf src/python/out
	rm -rf src/python/*.pyc
	$(MAKE) clean -C src/burkadt
	rm -rf $(BIN_INSTALL)/*
	rm -rf $(LIB_INSTALL)/*

