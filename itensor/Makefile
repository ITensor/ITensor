#
# Makefile for ITensor library
#
################################

include ../this_dir.mk
include ../options.mk

####################################

SOURCES=itensor.cc iqtensor.cc svdworker.cc mps.cc mpo.cc dmrg.cc

HEADERS=types.h allocator.h real.h permutation.h index.h prodstats.h \
        itensor.h iqindex.h iqtensor.h condenser.h combiner.h iqcombiner.h svdworker.h \
        mps.h mpo.h DMRGOpts.h dmrg.h core.h BaseDMRGOpts.h \
        BaseDMRGWorker.h DMRGWorker.h Sweeps.h hams.h measure.h model.h\
        hams/hubbardchain.h hams/heisenberg.h hams/triheisenberg.h\
        model/spinhalf.h model/spinone.h model/hubbard.h model/spinless.h\
        eigensolver.h projectedop.h

####################################

INCLUDEFLAGS=-I$(INCLUDEDIR) -I$(BOOST_DIR) 

CCFLAGS= $(INCLUDEFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= $(INCLUDEFLAGS) -ansi -Wall -pedantic -Wno-long-long -g -O0\
 -DDEBUG -DBOUNDS -DMATRIXBOUNDS -DSTRONG_DEBUG

OBJECTS= $(patsubst %.cc,%.o, $(SOURCES))
GOBJECTS= $(patsubst %.cc,.debug_objs/%.o, $(SOURCES))
INCLUDE_HEADERS= $(patsubst %,$(INCLUDEDIR)/%, $(HEADERS))

.SUFFIXES: .cc $(SUFFIXES)

$(INCLUDEDIR)/hams/%.h: hams/%.h
	mkdir -p $(INCLUDEDIR)/hams
	cp $^ $(INCLUDEDIR)/hams

$(INCLUDEDIR)/model/%.h: model/%.h
	mkdir -p $(INCLUDEDIR)/model
	cp $^ $(INCLUDEDIR)/model

$(INCLUDEDIR)/%.h: %.h
	cp $^ $(INCLUDEDIR)

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

headerscopy: $(INCLUDE_HEADERS)
	cp core.h $(INCLUDEDIR)

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

clean:	
	rm -fr *.o .debug_objs libitensor.a libitensor-g.a

DEPHEADERS=types.h allocator.h real.h permutation.h index.h itensor.h
itensor.o: $(DEPHEADERS)
.debug_objs/itensor.o: $(DEPHEADERS)
DEPHEADERS+= iqindex.h iqtensor.h
iqtensor.o: $(DEPHEADERS)
.debug_objs/iqtensor.o: $(DEPHEADERS)
DEPHEADERS+= combiner.h condenser.h iqcombiner.h svdworker.h
svdworker.o: $(DEPHEADERS)
DEPHEADERS+= mps.h
mps.o: $(DEPHEADERS)
.debug_objs/mps.o: $(DEPHEADERS)
DEPHEADERS+= mpo.h
mpo.o: $(DEPHEADERS)
.debug_objs/mpo.o: $(DEPHEADERS)
DEPHEADERS+= DMRGOpts.h dmrg.h
dmrg.o: $(DEPHEADERS)
.debug_objs/dmrg.o: $(DEPHEADERS)
