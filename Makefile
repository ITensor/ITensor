#
# Makefile for ITensor library
#
####################################

include options.mk

build: configure link_Makefiles
	cd utilities && make
	cd matrix/$(PLATFORM) && make
	cd itensor && make install

configure:
	@echo THIS_DIR=`pwd` |cat - options.mk > .tmp && mv .tmp options.mk

link_Makefiles:
	cd utilities && ln -f -s Makefile.default Makefile
	cd matrix/$(PLATFORM) && ln -f -s Makefile.default Makefile
	cd itensor && ln -f -s Makefile.default Makefile

clean: link_Makefiles
	cd utilities && make clean
	cd matrix/$(PLATFORM) && make clean
	cd itensor && make clean
	cd sample && make clean
	cd sandbox && make clean
	rm -f include/*
	rm -f lib/*

distclean: clean
	cd utilities && rm -f Makefile
	cd matrix/$(PLATFORM) && rm -f Makefile
	cd itensor && rm -f Makefile
	cd sample && rm -f Makefile
	cd sandbox && rm -f Makefile
