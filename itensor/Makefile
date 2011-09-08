#
# Makefile for ITensor library
#
################################

include ../this_dir.mk
include ../options.mk

####################################

OBJECTS=itensor.o iqtensor.o mps.o

SOURCES=itensor.cc iqtensor.cc mps.cc

####################################

INCLUDEFLAGS=-I$(INCLUDEDIR) -I$(BOOST_DIR)

CCFLAGS= $(INCLUDEFLAGS) -O2 -DNDEBUG -DBOOST_DISABLE_ASSERTS
CCGFLAGS= $(INCLUDEFLAGS) -ansi -Wall -pedantic -g -O0 -DBOUNDS -DITENSOR_USE_AT -Wno-long-long -DMATRIXBOUNDS

GOBJECTS= $(patsubst %,.debug_objs/%, $(OBJECTS))

.SUFFIXES: .cc $(SUFFIXES)

%.o: %.cc
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

install: installo installg

libitensor.a: $(OBJECTS)
	ar rv libitensor.a $(OBJECTS)
	ranlib libitensor.a

libitensor-g.a: mkdebugdir $(GOBJECTS)
	ar rv libitensor-g.a $(GOBJECTS)
	ranlib libitensor-g.a

headerscopy:	
	cp -r *.h .profiling $(INCLUDEDIR)

libs: libitensor.a 
	cp libitensor.a $(LIBDIR)

installo: headerscopy libitensor.a
	cp libitensor.a $(LIBDIR)
	ranlib $(LIBDIR)/libitensor.a

installg: headerscopy libitensor-g.a
	cp libitensor-g.a $(LIBDIR)
	ranlib $(LIBDIR)/libitensor-g.a

mkdebugdir:
	mkdir -p .debug_objs

itensor.o: types.h real.h permutation.h index.h itensor.h
iqtensor.o: types.h real.h permutation.h index.h itensor.h \
combiner.h iqindex.h iqcombiner.h iqtensor.h
mps.o: types.h real.h permutation.h index.h itensor.h \
combiner.h iqindex.h iqcombiner.h iqtensor.h model.h mps.h

.debug_objs/itensor.o: types.h real.h permutation.h index.h itensor.h
.debug_objs/iq.o: types.h real.h permutation.h index.h \
itensor.h combiner.h iqindex.h iqcombiner.h iqtensor.h 
.debug_objs/mps.o: types.h real.h permutation.h index.h \
itensor.h combiner.h iqindex.h iqcombiner.h iqtensor.h model.h mps.h

clean:	
	rm -fr *.o .debug_objs libitensor.a libitensor-g.a

#dependencies
