#
# Makefile for ITensor library
#
################################

include ../this_dir.mk
include ../options.mk

####################################

HEADERS=

OBJECTS=itensor.o

SOURCES=itensor.cc

PYTHON_INCLUDE=/opt/local/lib -I/opt/local/Library/Frameworks/Python.framework/Versions/2.6/include/python2.6
PYTHON_LIB=/opt/local/Library/Frameworks/Python.framework/Versions/2.6/lib
BOOST_LIBDIR=$(BOOST_DIR)/../lib

####################################

INCLUDEFLAGS=-I$(INCLUDEDIR) -I$(BOOST_DIR) -I$(PYTHON_INCLUDE)

CCFLAGS= $(INCLUDEFLAGS) -O2 -DNDEBUG -DBOOST_DISABLE_ASSERTS
CCGFLAGS= $(INCLUDEFLAGS) -ansi -g -O0 -DBOUNDS -DITENSOR_USE_AT

LIBDIRS=-L$(BOOST_LIBDIR) -L$(LIBDIR) -L$(PYTHON_LIB) 
LIBFLAGS=-lpython2.6 -lboost_python -litensor -lutilities -lmatrix -framework Accelerate -dynamiclib
LIBGFLAGS=-lpython2.6 -lboost_python -litensor-g -lutilities-g -lmatrix-g -framework Accelerate -dynamiclib

GOBJECTS= $(patsubst %,.debug_objs/%, $(OBJECTS))

.SUFFIXES: .cc $(SUFFIXES)

%.o: %.cc
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

build: itensor.o
	$(CCCOM) $(INCLUDEFLAGS) $(LIBDIRS) $(LIBFLAGS) itensor.o -o itensor.so

debug: mkdebugdir .debug_objs/itensor.o
	$(CCCOM) $(INCLUDEFLAGS) $(LIBDIRS) $(LIBGFLAGS) .debug_objs/itensor.o -o itensor.so

mkdebugdir:	
	mkdir -p .debug_objs

clean:	
	rm -fr *.o *.so .debug_objs

#dependencies
