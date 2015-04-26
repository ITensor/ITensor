#
# Makefile for ITensor library
#
################################

include ../this_dir.mk
include ../options.mk

####################################

SOURCES+= simpletensor.cc 
SOURCES+= matrix/vec.cc 
SOURCES+= matrix/mat.cc 
#SOURCES+= matrix/algs.cc 
#SOURCES = index.cc 
#SOURCES+= contract.cc 
#SOURCES+= itensor.cc 
#SOURCES+= iqindex.cc 
#SOURCES+= iqtensor.cc 

#SOURCES+= spectrum.cc 
#SOURCES+= svdalgs.cc 
#SOURCES+= mps.cc 
#SOURCES+= mpo.cc 
#SOURCES+= tevol.cc
#SOURCES+= autompo.cc

HEADERS=detail/algs.h detail/call_rewrite.h detail/printing.h detail/skip_iterator.h \
        detail/gcounter.h global.h real.h permutation.h index.h \
        lapack_wrap.h matrix/strideiter.h matrix/vec.h \
        matrix/mrange.h matrix/miterator.h matrix/mat.h matrix/algs.h \
        range.h indexset.h counter.h simpletensor.h contract.h \
        itdata/itdata.h \
        itdata/itreal.h itdata/itdiag.h itdata/itcombiner.h itdata/itcplx.h \
        itdata/iqtdata.h \
        itensor.h itensor.ih qn.h iqindex.h  \
        iqtensor.h iqtensor.ih \
        #condenser.h qcounter.h iqcombiner.h \
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

$(ITENSOR_INCLUDEDIR)/matrix/%.h: matrix/%.h
	mkdir -p $(ITENSOR_INCLUDEDIR)/matrix
	cp $^ $(ITENSOR_INCLUDEDIR)/matrix

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

$(ITENSOR_INCLUDEDIR)/%.ih: %.ih
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

debug: installg
installg: headerscopy libitensor-g.a
	cp libitensor-g.a $(ITENSOR_LIBDIR)
	ranlib $(ITENSOR_LIBDIR)/libitensor-g.a

mkdebugdir:
	mkdir -p .debug_objs
	mkdir -p .debug_objs/matrix

clean:	
	rm -fr *.o .debug_objs libitensor.a libitensor-g.a

GDEPHEADERS=real.h global.h index.h permutation.h
index.o: $(GDEPHEADERS)
.debug_objs/index.o: $(GDEPHEADERS)
GDEPHEADERS+= indexset.h
indexset.o: $(GDEPHEADERS)
.debug_objs/indexset.o: $(GDEPHEADERS)
GDEPHEADERS+= matrix/strideiter.h matrix/vec.h
matrix/vec.o: $(GDEPHEADERS)
.debug_objs/matrix/vec.o: $(GDEPHEADERS)
GDEPHEADERS+= matrix/mrange.h matrix/miterator.h matrix/mat.h
matrix/mat.o: $(GDEPHEADERS)
.debug_objs/matrix/mat.o: $(GDEPHEADERS)
GDEPHEADERS+= matrix/algs.h
matrix/algs.o: $(GDEPHEADERS)
.debug_objs/matrix/algs.o: $(GDEPHEADERS)
GDEPHEADERS+= contract.h range.h simpletensor.h
contract.o: $(GDEPHEADERS)
.debug_objs/contract.o: $(GDEPHEADERS)
ITDEPHEADERS= itensor.ih itensor.h itdata/itreal.h itdata/itdiag.h itdata/itcplx.h
GDEPHEADERS+= counter.h itdata/itdata.h
itensor.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itensor.o: $(ITDEPHEADERS) $(GDEPHEADERS)
IQDEPHEADERS= qn.h iqindex.h
iqindex.o: $(IQDEPHEADERS) $(GDEPHEADERS)
.debug_objs/iqindex.o: $(IQDEPHEADERS) $(GDEPHEADERS)
IQDEPHEADERS+= itdata/iqtdata.h iqtensor.ih iqtensor.h qcounter.h detail/skip_iterator.h
iqtensor.o: $(IQDEPHEADERS) $(GDEPHEADERS)
.debug_objs/iqtensor.o: $(IQDEPHEADERS) $(GDEPHEADERS)
