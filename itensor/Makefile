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
 -DDEBUG -DBOUNDS -DMATRIXBOUNDS -DSTRONG_DEBUG

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

HEADERS=types.h allocator.h real.h permutation.h index.h itensor.h
itensor.o: $(HEADERS)
.debug_objs/itensor.o: $(HEADERS)
HEADERS+= iqindex.h iqtensor.h
iqtensor.o: $(HEADERS)
.debug_objs/iqtensor.o: $(HEADERS)
HEADERS+= mps.h
mps.o: $(HEADERS)
.debug_objs/mps.o: $(HEADERS)
HEADERS+= mpo.h
mpo.o: $(HEADERS)
.debug_objs/mpo.o: $(HEADERS)
HEADERS+= dmrg.h
dmrg.o: $(HEADERS)
.debug_objs/dmrg.o: $(HEADERS)

clean:	
	rm -fr *.o .debug_objs libitensor.a libitensor-g.a

#dependencies
