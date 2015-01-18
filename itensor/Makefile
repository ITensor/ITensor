#
# Makefile for ITensor library
#
################################

include ../this_dir.mk
include ../options.mk

####################################

SOURCES = index.cc 
SOURCES+= itdata_functions.cc 
SOURCES+= simpletensor.cc 
SOURCES+= contract.cc 
SOURCES+= itensor.cc 
SOURCES+= iqindex.cc 
SOURCES+= iqtensor.cc 
#SOURCES+= condenser.cc
#SOURCES+= iqcombiner.cc 
#SOURCES+= spectrum.cc 
#SOURCES+= svdalgs.cc 
#SOURCES+= mps.cc 
#SOURCES+= mpo.cc 
#SOURCES+= tevol.cc
#SOURCES+= autompo.cc

HEADERS=detail/functions.h detail/skip_iterator.h detail/gcounter.h \
        global.h real.h permutation.h index.h \
        range.h indexset.h counter.h simpletensor.h contract.h \
        itdata/itdata.h itdata_functions.h \
        itdata/itdense.h itdata/itdiag.h itdata/itcombiner.h \
        itdata/iqtdata.h \
        itensor.h qn.h iqindex.h  \
        detail/skip_iterator.h iqtensor.h \
        condenser.h qcounter.h iqcombiner.h \
        spectrum.h svdalgs.h mps.h mpo.h core.h observer.h DMRGObserver.h \
        sweeps.h stats.h siteset.h\
        hams/HubbardChain.h hams/Heisenberg.h hams/ExtendedHubbard.h \
        hams/TriHeisenberg.h hams/Ising.h hams/J1J2Chain.h \
        hams/tJChain.h \
        sites/spinhalf.h sites/spinone.h sites/hubbard.h sites/spinless.h\
        sites/tj.h sites/Z3.h\
        eigensolver.h localop.h localmpo.h localmposet.h \
        partition.h hambuilder.h localmpo_mps.h tevol.h dmrg.h bondgate.h\
        integrators.h idmrg.h TEvolObserver.h iterpair.h autompo.h



####################################

CCFLAGS= $(ITENSOR_INCLUDEFLAGS) $(OPTIMIZATIONS) -DPLATFORM_$(PLATFORM)
CCGFLAGS= $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS) -DPLATFORM_$(PLATFORM)

OBJECTS= $(patsubst %.cc,%.o, $(SOURCES))
GOBJECTS= $(patsubst %.cc,.debug_objs/%.o, $(SOURCES))
INCLUDE_HEADERS= $(patsubst %,$(ITENSOR_INCLUDEDIR)/%, $(HEADERS))

.SUFFIXES: .cc $(SUFFIXES)

$(ITENSOR_INCLUDEDIR)/itdata/%.h: itdata/%.h
	mkdir -p $(ITENSOR_INCLUDEDIR)/itdata
	cp $^ $(ITENSOR_INCLUDEDIR)/itdata

$(ITENSOR_INCLUDEDIR)/detail/%.h: detail/%.h
	mkdir -p $(ITENSOR_INCLUDEDIR)/detail
	cp $^ $(ITENSOR_INCLUDEDIR)/detail

$(ITENSOR_INCLUDEDIR)/hams/%.h: hams/%.h
	mkdir -p $(ITENSOR_INCLUDEDIR)/hams
	cp $^ $(ITENSOR_INCLUDEDIR)/hams

$(ITENSOR_INCLUDEDIR)/sites/%.h: sites/%.h
	mkdir -p $(ITENSOR_INCLUDEDIR)/sites
	cp $^ $(ITENSOR_INCLUDEDIR)/sites

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

DEPHEADERS=real.h global.h index.h permutation.h
index.o: $(DEPHEADERS)
.debug_objs/index.o: $(DEPHEADERS)
DEPHEADERS+= indexset.h
indexset.o: $(DEPHEADERS)
.debug_objs/indexset.o: $(DEPHEADERS)
DEPHEADERS+= contract.h
contract.o: $(DEPHEADERS)
.debug_objs/contract.o: $(DEPHEADERS)
DEPHEADERS+= simpletensor.h itdata/itdata.h itdata_functions.h \
itensor.h counter.h
itensor.o: $(DEPHEADERS)
.debug_objs/itensor.o: $(DEPHEADERS)
DEPHEADERS+= qn.h iqindex.h
iqindex.o: $(DEPHEADERS)
.debug_objs/iqindex.o: $(DEPHEADERS)
DEPHEADERS+= itdata/iqtdata.h iqtensor.h qcounter.h detail/skip_iterator.h
iqtensor.o: $(DEPHEADERS)
.debug_objs/iqtensor.o: $(DEPHEADERS)
DEPHEADERS+= condenser.h iqcombiner.h localmpo.h
iqcombiner.o: $(DEPHEADERS)
.debug_objs/iqcombiner.o: $(DEPHEADERS)
condenser.o: $(DEPHEADERS)
.debug_objs/condenser.o: $(DEPHEADERS)
DEPHEADERS+= spectrum.h
spectrum.o: $(DEPHEADERS)
.debug_objs/spectrum.o: $(DEPHEADERS)
DEPHEADERS+= svdalgs.h
svdalgs.o: $(DEPHEADERS)
.debug_objs/svdalgs.o: $(DEPHEADERS)
DEPHEADERS+= mps.h
mps.o: $(DEPHEADERS)
.debug_objs/mps.o: $(DEPHEADERS)
DEPHEADERS+= mpo.h
mpo.o: $(DEPHEADERS)
.debug_objs/mpo.o: $(DEPHEADERS)
DEPHEADERS+= autompo.h
autompo.o: $(DEPHEADERS)
.debug_objs/autompo.o: $(DEPHEADERS)
DEPHEADERS+= bondgate.h tevol.h integrators.h
tevol.o: $(DEPHEADERS)
.debug_objs/tevol.o: $(DEPHEADERS)
