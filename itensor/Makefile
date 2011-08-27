#
# Makefile for ITensor library
#
################################

include ../this_dir.mk
include ../options.mk

####################################

HEADERS=tensor.h mps.h algorithm.h iq.h iqcombiner.h hams.h

OBJECTS=tensor.o mps.o algorithm.o iq.o

SOURCES=tensor.cc mps.cc algorithm.cc iq.cc

####################################

INCLUDEFLAGS=-I$(INCLUDEDIR) -I$(BOOST_DIR)

CCFLAGS= $(INCLUDEFLAGS) -O2 -DNDEBUG -DBOOST_DISABLE_ASSERTS
CCGFLAGS= $(INCLUDEFLAGS) -ansi -g -O0 -DBOUNDS -DITENSOR_USE_AT

GOBJECTS= $(patsubst %,.debug_objs/%, $(OBJECTS))

.SUFFIXES: .cc $(SUFFIXES)

%.o: %.cc
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

install: installo installg

libitensor.a: $(OBJECTS) $(HEADERS)
	ar rv libitensor.a $(OBJECTS)
	ranlib libitensor.a

libitensor-g.a: mkdebugdir $(GOBJECTS) $(HEADERS)
	ar rv libitensor-g.a $(GOBJECTS)
	ranlib libitensor-g.a

headerscopy:	
	cp $(HEADERS) $(INCLUDEDIR)

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

tensor.o: tensor.h
iq.o: tensor.h iqcombiner.h iq.h
mps.o: tensor.h iqcombiner.h iq.h mps.h
algorithm.o: tensor.h iqcombiner.h iq.h mps.h algorithm.h 

.debug_objs/tensor.o: tensor.h
.debug_objs/iq.o: tensor.h iqcombiner.h iq.h
.debug_objs/mps.o: tensor.h iqcombiner.h iq.h mps.h
.debug_objs/algorithm.o: tensor.h iqcombiner.h iq.h mps.h algorithm.h 

clean:	
	rm -fr *.o .debug_objs libitensor.a libitensor-g.a

#dependencies
