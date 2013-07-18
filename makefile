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

binary = $(BIN_INSTALL)/$(progname)
sourcecode = $(shell ls src/ctqw/*.f90)
LibmatExp = $(LIB_INSTALL)/libr8$(libName).a $(LIB_INSTALL)/libc8$(libName).a $(LIB_INSTALL)/libmatExp$(libName).a

#Makefile
all: $(binary) python

$(binary): $(LibmatExp) $(sourcecode)
	cd src/ctqw; sed 's/^\tuse/!\tuse/g;s/dbesjn/bessel_jn/g;/^\tsubroutine progressbar.*/,/^\tend subroutine progressbar.*/s/^/!/;/call progressbar(m,terms)/s/^/!/' libctqw_intel.f90 > libctqw_gcc.f90
	$(MAKE) -C src/ctqw
	
$(LibmatExp):
	$(MAKE) -C src/burkadt

$(sourcecode):
	
# Python wrappers
python: src/ctqw/libctqw$(libName).f90 $(LibmatExp) 
	$(MAKE) python -C src/ctqw

#Cleaning files
clean:
	$(MAKE) clean -C src/ctqw
	$(MAKE) clean -C src/burkadt
	rm -f $(binary)

fullclean: clean
	rm -rf $(BIN_INSTALL)/*
	rm -rf $(LIB_INSTALL)/*

