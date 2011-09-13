#
# Makefile for ITensor library
#
################################

include ../this_dir.mk
include ../options.mk

####################################

SOURCES=itensor.cc iqtensor.cc mps.cc mpo.cc dmrg.cc

####################################

INCLUDEFLAGS=-I$(INCLUDEDIR) -I$(BOOST_DIR)

CCFLAGS= $(INCLUDEFLAGS) -O2
CCGFLAGS= $(INCLUDEFLAGS) -ansi -Wall -pedantic -Wno-long-long -g -O0\
 -DDEBUG -DBOUNDS -DMATRIXBOUNDS

OBJECTS= $(patsubst %.cc,%.o, $(SOURCES))
GOBJECTS= $(patsubst %.cc,.debug_objs/%.o, $(SOURCES))

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

itensor.o: types.h allocator.h real.h permutation.h index.h itensor.h
iqtensor.o: itensor.o iqindex.h iqtensor.h
mps.o: iqtensor.o mps.h
mpo.o: mps.o mpo.h
dmrg.o: mpo.o dmrg.h

.debug_objs/itensor.o: types.h allocator.h real.h permutation.h index.h itensor.h
.debug_objs/iqtensor.o: .debug_objs/itensor.o iqindex.h iqtensor.h
.debug_objs/mps.o: .debug_objs/iqtensor.o mps.h
.debug_objs/mpo.o: .debug_objs/mps.o mpo.h
.debug_objs/dmrg.o: .debug_objs/mpo.o dmrg.h

clean:	
	rm -fr *.o .debug_objs libitensor.a libitensor-g.a

#dependencies
