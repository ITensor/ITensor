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
	@echo "THIS_DIR=$(PWD)" > this_dir.mk
	@echo "#ifndef __ITENSOR_CONFIG_H" > itensor/config.h
	@echo "#define __ITENSOR_CONFIG_H" >> itensor/config.h
	@echo "" >> itensor/config.h
	@echo "#ifndef PLATFORM_$(PLATFORM)" >> itensor/config.h
	@echo "#define PLATFORM_$(PLATFORM)" >> itensor/config.h
	@echo "#endif" >> itensor/config.h
	@echo "" >> itensor/config.h
	@echo "#ifndef __ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES" >> itensor/config.h
	@echo "#define __ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES 0" >> itensor/config.h
	@echo "#endif" >> itensor/config.h
	@echo "" >> itensor/config.h
ifdef ITENSOR_USE_HDF5
ifeq ($(ITENSOR_USE_HDF5),1)
	@echo "#ifndef ITENSOR_USE_HDF5" >> itensor/config.h
	@echo "#define ITENSOR_USE_HDF5 1" >> itensor/config.h
	@echo "#endif" >> itensor/config.h
endif
endif
	@echo "#endif " >> itensor/config.h

clean:
	@echo "Removing temporary build files"
	@touch this_dir.mk
	@cd itensor && $(MAKE) clean
	@cd sample && $(MAKE) clean
	@cd unittest && $(MAKE) clean
	@rm -f lib/*
	@rm -f this_dir.mk
	@rm -f itensor/config.h

distclean: clean
	@rm -f this_dir.mk options.mk
