#
# Makefile for ITensor library
#
####################################

include options.mk

build: configure 
	cd utilities && make
	cd matrix/$(PLATFORM) && make
	cd itensor && make

configure:
	@echo THIS_DIR=`pwd` > this_dir.mk

clean:
	cd utilities && make clean
	cd matrix/$(PLATFORM) && make clean
	cd itensor && make clean
	cd sample && make clean
	cd sandbox && make clean
	rm -f include/*
	rm -f lib/*
