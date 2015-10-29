#
# Makefile for ITensor library
#
################################

include ../this_dir.mk
include ../options.mk

####################################

UTILSOURCES+= util/error.cc 
UTILSOURCES+= util/ran1.cc 
UTILSOURCES+= util/tarray1.cc 
UTILSOURCES+= util/input.cc 
UTILSOURCES+= util/option.cc 
UTILSOURCES+= util/cputime.cc 

MATSOURCES+= matrix/matrix.cc 
MATSOURCES+= matrix/utility.cc 
MATSOURCES+= matrix/sparse.cc 
MATSOURCES+= matrix/david.cc 
MATSOURCES+= matrix/matrixref.cc 
MATSOURCES+= matrix/storelink.cc 
MATSOURCES+= matrix/sparseref.cc 
MATSOURCES+= matrix/daxpy.cc 
MATSOURCES+= matrix/svd.cc 

SOURCES+= index.cc 
SOURCES+= itensor.cc 
SOURCES+= iqindex.cc 
SOURCES+= iqtensor.cc 
SOURCES+= condenser.cc
SOURCES+= iqcombiner.cc 
SOURCES+= spectrum.cc 
SOURCES+= svdalgs.cc 
SOURCES+= mps/mps.cc 
SOURCES+= mps/mpo.cc 
SOURCES+= mps/tevol.cc
SOURCES+= mps/autompo.cc

HEADERS=util/tinyformat.h util/print.h util/types.h util/indent.h util/cputime.h util/tarray1.h \
        util/flstring.h util/option.h util/input.h util/error.h util/minmax.h util/cplx_literal.h \
        util/autovector.h util/optional_ptr.h util/count.h util/stdx.h util/timers.h \
        matrix/matrixref.h matrix/matrix.h matrix/sparse.h matrix/bigmatrix.h matrix/davidson.h\
	    matrix/storelink.h matrix/conjugate_gradient.h matrix/sparseref.h\
        matrix/svd.h \
        global.h real.h permutation.h indextype.h index.h \
        indexset.h counter.h simplematrix.h itensor.h itensor.ih \
        qn.h iqindex.h iqtdat.h  \
        detail/skip_iterator.h iqtensor.h iqtensor.ih \
        condenser.h combiner.h qcounter.h iqcombiner.h \
        spectrum.h svdalgs.h mps/mps.h mps/mpo.h core.h mps/observer.h mps/DMRGObserver.h \
        mps/sweeps.h util/stats.h mps/siteset.h\
        mps/hams/HubbardChain.h mps/hams/Heisenberg.h mps/hams/ExtendedHubbard.h \
        mps/hams/TriHeisenberg.h mps/hams/Ising.h mps/hams/J1J2Chain.h \
        mps/hams/tJChain.h \
        mps/sites/spinhalf.h mps/sites/spinone.h mps/sites/hubbard.h mps/sites/spinless.h\
        mps/sites/tj.h mps/sites/Z3.h\
        eigensolver.h mps/localop.h mps/localmpo.h mps/localmposet.h \
        partition.h mps/hambuilder.h mps/localmpo_mps.h mps/tevol.h mps/dmrg.h mps/bondgate.h\
        mps/integrators.h mps/idmrg.h mps/TEvolObserver.h iterpair.h mps/autompo.h \
        mps/lattice/latticebond.h mps/lattice/square.h mps/lattice/triangular.h



####################################

CCFLAGS= $(ITENSOR_INCLUDEFLAGS) $(OPTIMIZATIONS) -DPLATFORM_$(PLATFORM)
CCGFLAGS= $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS) -DPLATFORM_$(PLATFORM)

OBJECTS= $(patsubst %.cc,%.o, $(SOURCES))
GOBJECTS= $(patsubst %.cc,.debug_objs/%.o, $(SOURCES))
UTILOBJECTS= $(patsubst %.cc,%.o, $(UTILSOURCES))
UTILGOBJECTS= $(patsubst %.cc,.debug_objs/%.o, $(UTILSOURCES))
MATOBJECTS= $(patsubst %.cc,%.o, $(MATSOURCES))
MATGOBJECTS= $(patsubst %.cc,.debug_objs/%.o, $(MATSOURCES))
INCLUDE_HEADERS= $(patsubst %,$(ITENSOR_INCLUDEDIR)/%, $(HEADERS))

.SUFFIXES: .cc $(SUFFIXES)

$(ITENSOR_INCLUDEDIR)/util/%.h: util/%.h
	@cp $^ $(ITENSOR_INCLUDEDIR)

$(ITENSOR_INCLUDEDIR)/matrix/%.h: matrix/%.h
	@cp $^ $(ITENSOR_INCLUDEDIR)

$(ITENSOR_INCLUDEDIR)/detail/%.h: detail/%.h
	@mkdir -p $(ITENSOR_INCLUDEDIR)/detail
	@cp $^ $(ITENSOR_INCLUDEDIR)/detail

$(ITENSOR_INCLUDEDIR)/mps/hams/%.h: mps/hams/%.h
	@mkdir -p $(ITENSOR_INCLUDEDIR)/hams
	@cp $^ $(ITENSOR_INCLUDEDIR)/hams

$(ITENSOR_INCLUDEDIR)/mps/lattice/%.h: mps/lattice/%.h
	@mkdir -p $(ITENSOR_INCLUDEDIR)/lattice
	@cp $^ $(ITENSOR_INCLUDEDIR)/lattice

$(ITENSOR_INCLUDEDIR)/mps/sites/%.h: mps/sites/%.h
	@mkdir -p $(ITENSOR_INCLUDEDIR)/sites
	@cp $^ $(ITENSOR_INCLUDEDIR)/sites

$(ITENSOR_INCLUDEDIR)/mps/%.h: %.h
	@cp $^ $(ITENSOR_INCLUDEDIR)

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

libutilities.a:	$(UTILOBJECTS) $(HEADERS)
	@echo "Building static library libutilities.a"
	@ar r libutilities.a $(UTILOBJECTS)
	@ranlib libutilities.a

libutilities-g.a: mkdebugdir $(UTILGOBJECTS) $(HEADERS)
	@echo "Building static library libutilities-g.a"
	@ar r libutilities-g.a $(UTILGOBJECTS)
	@ranlib libutilities-g.a

libmatrix.a: $(MATOBJECTS) $(HEADERS)
	@echo "Building static library libmatrix.a"
	@ar r libmatrix.a $(MATOBJECTS)
	@touch libmatrix.a

libmatrix-g.a: mkdebugdir $(MATGOBJECTS) $(HEADERS)
	@echo "Building static library libmatrix-g.a"
	@ar r libmatrix-g.a $(MATGOBJECTS)
	@touch libmatrix-g.a

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

installo: headerscopy libutilities.a libmatrix.a libitensor.a
	@echo "Copying libutilities.a to $(ITENSOR_LIBDIR)"
	@cp libutilities.a $(ITENSOR_LIBDIR)
	@echo "Copying libmatrix.a to $(ITENSOR_LIBDIR)"
	@cp libmatrix.a $(ITENSOR_LIBDIR)
	@echo "Copying libitensor.a to $(ITENSOR_LIBDIR)"
	@cp libitensor.a $(ITENSOR_LIBDIR)
	@ranlib $(ITENSOR_LIBDIR)/libitensor.a

installg: headerscopy libutilities-g.a libmatrix-g.a libitensor-g.a
	@echo "Copying libutilities-g.a to $(ITENSOR_LIBDIR)"
	@cp libutilities-g.a $(ITENSOR_LIBDIR)
	@echo "Copying libmatrix-g.a to $(ITENSOR_LIBDIR)"
	@cp libmatrix-g.a $(ITENSOR_LIBDIR)
	@echo "Copying libitensor-g.a to $(ITENSOR_LIBDIR)"
	@cp libitensor-g.a $(ITENSOR_LIBDIR)
	@ranlib $(ITENSOR_LIBDIR)/libitensor-g.a

mkdebugdir:
	@mkdir -p .debug_objs
	@mkdir -p .debug_objs/util
	@mkdir -p .debug_objs/matrix
	@mkdir -p .debug_objs/mps

clean:	
	@echo "Cleaning build files"
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
DEPHEADERS+= combiner.h condenser.h iqcombiner.h mps/localmpo.h
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
DEPHEADERS+= mps/mps.h
mps/mps.o: $(DEPHEADERS)
.debug_objs/mps.o: $(DEPHEADERS)
DEPHEADERS+= mps/mpo.h
mps/mpo.o: $(DEPHEADERS)
.debug_objs/mpo.o: $(DEPHEADERS)
DEPHEADERS+= mps/autompo.h
mps/autompo.o: $(DEPHEADERS)
.debug_objs/autompo.o: $(DEPHEADERS)
DEPHEADERS+= mps/bondgate.h mps/tevol.h mps/integrators.h
mps/tevol.o: $(DEPHEADERS)
.debug_objs/tevol.o: $(DEPHEADERS)
