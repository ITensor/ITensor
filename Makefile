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
	#@echo THIS_DIR=`pwd` > .tmp1
	#@tail -n +2 options.mk > .tmp2
	#@cat .tmp1 .tmp2 > options.mk
	#@rm .tmp1 .tmp2

clean:
	cd utilities && make clean
	cd matrix/$(PLATFORM) && make clean
	cd itensor && make clean
	cd sample && make clean
	cd sandbox && make clean
	rm -f include/*
	rm -f lib/*
