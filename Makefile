#
# Makefile for ITensor libraries
#
####################################

include options.mk

build: itensor 

itensor: configure
	@echo
	@echo Building ITensor library
	@echo
	@cd itensor && $(MAKE)

configure:
	@echo
	@echo Configure: Writing current dir to this_dir.mk
	@echo THIS_DIR=`pwd` > this_dir.mk

clean:
	@echo "Removing temporary build files"
	@touch this_dir.mk
	@cd itensor && $(MAKE) clean
	@cd sample && $(MAKE) clean
	@cd unittest && $(MAKE) clean
	@rm -f lib/*
	@rm -f this_dir.mk

distclean: clean
	@rm -f this_dir.mk options.mk
