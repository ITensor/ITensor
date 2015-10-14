#
# Makefile for ITensor library
#
################################

include ../this_dir.mk
include ../options.mk

####################################

SOURCES+= util/error.cc 
SOURCES+= util/ran1.cc 
SOURCES+= util/tarray1.cc 
SOURCES+= util/input.cc 
SOURCES+= util/option.cc 
SOURCES+= util/cputime.cc 

#SOURCES+= matrix/matrix.cc 
#SOURCES+= matrix/utility.cc 
#SOURCES+= matrix/sparse.cc 
#SOURCES+= matrix/david.cc 
#SOURCES+= matrix/matrixref.cc 
#SOURCES+= matrix/storelink.cc 
#SOURCES+= matrix/sparseref.cc 
#SOURCES+= matrix/daxpy.cc 
#SOURCES+= matrix/svd.cc 

#SOURCES+= index.cc 
#SOURCES+= itensor.cc 
#SOURCES+= iqindex.cc 
#SOURCES+= iqtensor.cc 
#SOURCES+= condenser.cc
#SOURCES+= iqcombiner.cc 
#SOURCES+= spectrum.cc 
#SOURCES+= svdalgs.cc 
#SOURCES+= mps.cc 
#SOURCES+= mpo.cc 
#SOURCES+= tevol.cc
#SOURCES+= autompo.cc

HEADERS= util/tinyformat.h util/print.h util/types.h util/indent.h util/cputime.h util/tarray1.h \
        util/flstring.h util/option.h util/input.h util/error.h util/minmax.h util/cplx_literal.h \
        util/autovector.h util/optional_ptr.h util/count.h util/stdx.h util/timers.h \
        global.h real.h permutation.h indextype.h index.h \
        indexset.h counter.h simplematrix.h itensor.h itensor.ih \
        qn.h iqindex.h iqtdat.h  \
        detail/skip_iterator.h iqtensor.h iqtensor.ih \
        condenser.h combiner.h qcounter.h iqcombiner.h \
        spectrum.h svdalgs.h mps.h mpo.h core.h observer.h DMRGObserver.h \
        sweeps.h stats.h siteset.h\
        hams/HubbardChain.h hams/Heisenberg.h hams/ExtendedHubbard.h \
        hams/TriHeisenberg.h hams/Ising.h hams/J1J2Chain.h \
        hams/tJChain.h \
        sites/spinhalf.h sites/spinone.h sites/hubbard.h sites/spinless.h\
        sites/tj.h sites/Z3.h\
        eigensolver.h localop.h localmpo.h localmposet.h \
        partition.h hambuilder.h localmpo_mps.h tevol.h dmrg.h bondgate.h\
        integrators.h idmrg.h TEvolObserver.h iterpair.h autompo.h \
        lattice/latticebond.h lattice/square.h lattice/triangular.h



####################################

CCFLAGS= $(ITENSOR_INCLUDEFLAGS) $(OPTIMIZATIONS) -DPLATFORM_$(PLATFORM)
CCGFLAGS= $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS) -DPLATFORM_$(PLATFORM)

OBJECTS= $(patsubst %.cc,%.o, $(SOURCES))
GOBJECTS= $(patsubst %.cc,.debug_objs/%.o, $(SOURCES))
INCLUDE_HEADERS= $(patsubst %,$(ITENSOR_INCLUDEDIR)/%, $(HEADERS))

.SUFFIXES: .cc $(SUFFIXES)

$(ITENSOR_INCLUDEDIR)/util/%.h: util/%.h
	@cp $^ $(ITENSOR_INCLUDEDIR)

$(ITENSOR_INCLUDEDIR)/matrix/%.h: matrix/%.h
	@cp $^ $(ITENSOR_INCLUDEDIR)

$(ITENSOR_INCLUDEDIR)/detail/%.h: detail/%.h
	@mkdir -p $(ITENSOR_INCLUDEDIR)/detail
	@cp $^ $(ITENSOR_INCLUDEDIR)/detail

$(ITENSOR_INCLUDEDIR)/hams/%.h: hams/%.h
	@mkdir -p $(ITENSOR_INCLUDEDIR)/hams
	@cp $^ $(ITENSOR_INCLUDEDIR)/hams

$(ITENSOR_INCLUDEDIR)/lattice/%.h: lattice/%.h
	@mkdir -p $(ITENSOR_INCLUDEDIR)/lattice
	@cp $^ $(ITENSOR_INCLUDEDIR)/lattice

$(ITENSOR_INCLUDEDIR)/sites/%.h: sites/%.h
	@mkdir -p $(ITENSOR_INCLUDEDIR)/sites
	@cp $^ $(ITENSOR_INCLUDEDIR)/sites

$(ITENSOR_INCLUDEDIR)/%.h: %.h
	@cp $^ $(ITENSOR_INCLUDEDIR)

$(ITENSOR_INCLUDEDIR)/%.ih: %.ih
	@cp $^ $(ITENSOR_INCLUDEDIR)

%.o: %.cc
	@echo "Compiling itensor/$< with optimizations"
	$(eval COMMAND = $(CCCOM) -c $(CCFLAGS) -o $@ $<)
	@$(COMMAND) || (echo "Failure while executing command: $(COMMAND)" && exit 1)

.debug_objs/%.o: %.cc
	@echo "Compiling itensor/$< in debug mode"
	$(eval COMMAND = $(CCCOM) -c $(CCGFLAGS) -o $@ $<)
	@$(COMMAND) || (echo "Failure while executing command: $(COMMAND)" && exit 1)

install: installo installg

libitensor.a: $(OBJECTS)
	@echo "Building static library libitensor.a"
	@ar r libitensor.a $(OBJECTS)
	@ranlib libitensor.a

libitensor-g.a: mkdebugdir $(GOBJECTS)
	@echo "Building static library libitensor-g.a"
	@ar r libitensor-g.a $(GOBJECTS)
	@ranlib libitensor-g.a

headerscopy: $(INCLUDE_HEADERS)
	@echo "Copying header files to $(ITENSOR_INCLUDEDIR)"
	@cp core.h $(ITENSOR_INCLUDEDIR)

libs: libitensor.a 
	@echo "Copying libitensor.a to $(ITENSOR_LIBDIR)"
	@cp libitensor.a $(ITENSOR_LIBDIR)

installo: headerscopy libitensor.a
	@echo "Copying libitensor.a to $(ITENSOR_LIBDIR)"
	@cp libitensor.a $(ITENSOR_LIBDIR)
	@ranlib $(ITENSOR_LIBDIR)/libitensor.a

installg: headerscopy libitensor-g.a
	@echo "Copying libitensor-g.a to $(ITENSOR_LIBDIR)"
	@cp libitensor-g.a $(ITENSOR_LIBDIR)
	@ranlib $(ITENSOR_LIBDIR)/libitensor-g.a

mkdebugdir:
	@mkdir -p .debug_objs
	@mkdir -p .debug_objs/util
	@mkdir -p .debug_objs/matrix
	@mkdir -p .debug_objs/mps

clean:	
	@echo "Cleaning build files in folder itensor/"
	@rm -fr *.o .debug_objs util/*.o matrix/*.o mps/*.o libitensor.a libitensor-g.a

DEPHEADERS=real.h global.h indextype.h index.h permutation.h
index.o: $(DEPHEADERS)
.debug_objs/index.o: $(DEPHEADERS)
DEPHEADERS+= indexset.h
indexset.o: $(DEPHEADERS)
.debug_objs/indexset.o: $(DEPHEADERS)
DEPHEADERS+=  itensor.h counter.h simplematrix.h
itensor.o: $(DEPHEADERS)
.debug_objs/itensor.o: $(DEPHEADERS)
DEPHEADERS+= qn.h iqindex.h
iqindex.o: $(DEPHEADERS)
.debug_objs/iqindex.o: $(DEPHEADERS)
DEPHEADERS+= iqtdat.h iqtensor.h qcounter.h detail/skip_iterator.h
iqtensor.o: $(DEPHEADERS)
.debug_objs/iqtensor.o: $(DEPHEADERS)
DEPHEADERS+= combiner.h condenser.h iqcombiner.h localmpo.h
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
