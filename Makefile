#
# Makefile for ITensor libraries
#
####################################

.NOTPARALLEL:

include options.mk

build: itensor 

itensor: configure
	@echo
	@echo Building ITensor library
	@echo
	cd itensor && make

configure:
	@echo
	@echo Configure: Writing current dir to this_dir.mk
	@echo THIS_DIR=`pwd` > this_dir.mk

clean:
	cd itensor && make clean
	cd sample && make clean
	cd unittest && make clean
	rm -f lib/*
	rm -f this_dir.mk

distclean: clean
	rm -f this_dir.mk options.mk
