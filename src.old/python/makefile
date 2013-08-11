#!/usr/bin/make -f

BIN_INSTALL = bin
BIN_PATH = ../../$(BIN_INSTALL)
LIB_INSTALL = lib
LIB_PATH = ../../$(LIB_INSTALL)

PROG = $(BIN_PATH)/pyctqw
CONFIG = $(BIN_PATH)/input.pyqw

#Makefile
all: $(PROG) $(CONFIG)

$(PROG): pyctqw
	cp pyctqw $(BIN_PATH)
	
$(CONFIG): input.pyqw
	cp input.pyqw $(BIN_PATH)

#Cleaning files
clean:
	rm -f *.pyc
	rm -f $(PROG) $(CONFIG)













