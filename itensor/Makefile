#
# Makefile for ITensor library
#
################################

include ../this_dir.mk
include ../options.mk

####################################

SOURCES=index.cc itensor.cc itsparse.cc \
        iqindex.cc iqtensor.cc iqcombiner.cc iqtsparse.cc\
        svdalgs.cc mps.cc mpo.cc tevol.cc

HEADERS=global.h allocator.h real.h permutation.h index.h prodstats.h \
        indexset.h counter.h itensor.h qn.h iqindex.h iqtensor.h \
        condenser.h combiner.h qcounter.h iqcombiner.h \
        spectrum.h svdalgs.h mps.h mpo.h core.h observer.h DMRGObserver.h \
        sweeps.h stats.h model.h\
        hams/HubbardChain.h hams/Heisenberg.h hams/ExtendedHubbard.h \
        hams/TriHeisenberg.h hams/Ising.h hams/J1J2Chain.h \
        hams/tJChain.h \
        model/spinhalf.h model/spinone.h model/hubbard.h model/spinless.h\
        model/tj.h \
        eigensolver.h localop.h localmpo.h localmposet.h itsparse.h iqtsparse.h\
        partition.h option.h hambuilder.h localmpo_mps.h tevol.h dmrg.h bondgate.h

####################################

CCFLAGS= $(ITENSOR_INCLUDEFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS)

OBJECTS= $(patsubst %.cc,%.o, $(SOURCES))
GOBJECTS= $(patsubst %.cc,.debug_objs/%.o, $(SOURCES))
INCLUDE_HEADERS= $(patsubst %,$(ITENSOR_INCLUDEDIR)/%, $(HEADERS))

.SUFFIXES: .cc $(SUFFIXES)

$(ITENSOR_INCLUDEDIR)/hams/%.h: hams/%.h
	mkdir -p $(ITENSOR_INCLUDEDIR)/hams
	cp $^ $(ITENSOR_INCLUDEDIR)/hams

$(ITENSOR_INCLUDEDIR)/model/%.h: model/%.h
	mkdir -p $(ITENSOR_INCLUDEDIR)/model
	cp $^ $(ITENSOR_INCLUDEDIR)/model

$(ITENSOR_INCLUDEDIR)/%.h: %.h
	cp $^ $(ITENSOR_INCLUDEDIR)

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
	cp core.h $(ITENSOR_INCLUDEDIR)

libs: libitensor.a 
	cp libitensor.a $(ITENSOR_LIBDIR)

installo: headerscopy libitensor.a
	cp libitensor.a $(ITENSOR_LIBDIR)
	ranlib $(ITENSOR_LIBDIR)/libitensor.a

installg: headerscopy libitensor-g.a
	cp libitensor-g.a $(ITENSOR_LIBDIR)
	ranlib $(ITENSOR_LIBDIR)/libitensor-g.a

mkdebugdir:
	mkdir -p .debug_objs

clean:	
	rm -fr *.o .debug_objs libitensor.a libitensor-g.a

DEPHEADERS=global.h real.h permutation.h index.h option.h
index.o: $(DEPHEADERS)
.debug_objs/index.o: $(DEPHEADERS)
DEPHEADERS+= indexset.h
indexset.o: $(DEPHEADERS)
.debug_objs/indexset.o: $(DEPHEADERS)
DEPHEADERS+= allocator.h itensor.h counter.h
itensor.o: $(DEPHEADERS)
.debug_objs/itensor.o: $(DEPHEADERS)
DEPHEADERS+= itsparse.h
itsparse.o: $(DEPHEADERS)
.debug_objs/itsparse.o: $(DEPHEADERS)
DEPHEADERS+= qn.h iqindex.h
iqindex.o: $(DEPHEADERS)
.debug_objs/iqindex.o: $(DEPHEADERS)
DEPHEADERS+= iqtensor.h qcounter.h
iqtensor.o: $(DEPHEADERS)
.debug_objs/iqtensor.o: $(DEPHEADERS)
DEPHEADERS+= iqtsparse.h
iqtsparse.o: $(DEPHEADERS)
.debug_objs/iqtsparse.o: $(DEPHEADERS)
DEPHEADERS+= combiner.h condenser.h iqcombiner.h localmpo.h
iqcombiner.o: $(DEPHEADERS)
.debug_objs/iqcombiner.o: $(DEPHEADERS)
DEPHEADERS+= spectrum.h svdalgs.h
svdalgs.o: $(DEPHEADERS)
.debug_objs/svdalgs.o: $(DEPHEADERS)
DEPHEADERS+= mps.h
mps.o: $(DEPHEADERS)
.debug_objs/mps.o: $(DEPHEADERS)
DEPHEADERS+= mpo.h
mpo.o: $(DEPHEADERS)
.debug_objs/mpo.o: $(DEPHEADERS)
DEPHEADERS+= bondgate.h tevol.h integrators.h
tevol.o: $(DEPHEADERS)
.debug_objs/tevol.o: $(DEPHEADERS)
