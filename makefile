#!/usr/bin/make -f

progname = ctqw
BIN_INSTALL = bin
LIB_INSTALL = lib

ifeq ($(FC),gfortran)
	libName = _gcc
else
	FC = ifort
	libName = _intel
endif

binary = $(BIN_INSTALL)/$(progname)
LibmatExp = $(LIB_INSTALL)/libr8$(libName).a $(LIB_INSTALL)/libc8$(libName).a $(LIB_INSTALL)/libmatExp$(libName).a

#Makefile
$(binary): $(LibmatExp)
	$(MAKE) -C src/ctqw
	
$(LibmatExp):
	$(MAKE) -C src/burkadt

#Cleaning files
clean:
	$(MAKE) clean -C src/ctqw
	$(MAKE) clean -C src/burkadt
	rm $(binary)

fullclean:
	$(MAKE) clean -C src/ctqw
	$(MAKE) clean -C src/burkadt
	rm -rf $(BIN_INSTALL)/*
	rm -rf $(LIB_INSTALL)/*

