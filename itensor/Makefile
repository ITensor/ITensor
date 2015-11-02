#
# Makefile for ITensor library
#

include ../this_dir.mk
include ../options.mk

####################################

SOURCES+= util/args.cc     
SOURCES+= util/input.cc
SOURCES+= util/cputime.cc
SOURCES+= tensor/lapack_wrap.cc 
SOURCES+= tensor/vec.cc 
SOURCES+= tensor/mat.cc 
SOURCES+= tensor/gemm.cc 
SOURCES+= tensor/algs.cc 
SOURCES+= tensor/contract.cc 
SOURCES+= itdata/dense.cc 
SOURCES+= itdata/itcombiner.cc 
SOURCES+= itdata/diag.cc 
#SOURCES+= itdata/iqtreal.cc
#SOURCES+= itdata/iqtcombiner.cc 
#SOURCES+= itdata/iqtdiag.cc
##SOURCES+= itdata/itlazy.cc
SOURCES+= index.cc 
SOURCES+= itensor.cc 
#SOURCES+= iqindex.cc 
#SOURCES+= iqtensor.cc 
#SOURCES+= spectrum.cc 
#SOURCES+= svdalgs.cc 
#SOURCES+= mps/mps.cc 
#SOURCES+= mps/mpo.cc 
#SOURCES+= mps/autompo.cc

####################################

CCFLAGS= $(ITENSOR_INCLUDEFLAGS) $(OPTIMIZATIONS) -DPLATFORM_$(PLATFORM) -D__ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES=0
CCGFLAGS= $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS) -DPLATFORM_$(PLATFORM) -D__ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES=0

OBJECTS= $(patsubst %.cc,%.o, $(SOURCES))
GOBJECTS= $(patsubst %.cc,.debug_objs/%.o, $(SOURCES))

.SUFFIXES: .cc $(SUFFIXES)

%.o: %.cc
	@echo "Compiling itensor/$< with optimizations"
	$(eval COMMAND = $(CCCOM) -c $(CCFLAGS) -o $@ $<)
	@$(COMMAND) || (echo "Failure while executing command: $(COMMAND)" && exit 1)

.debug_objs/%.o: %.cc
	@echo "Compiling itensor/$< in debug mode"
	$(eval COMMAND = $(CCCOM) -c $(CCGFLAGS) -o $@ $<)
	@$(COMMAND) || (echo "Failure while executing command: $(COMMAND)" && exit 1)

install: build debug

libitensor.a: $(OBJECTS)
	@echo "Building static library $(ITENSOR_LIBDIR)/libitensor.a"
	@ar r $(ITENSOR_LIBDIR)/libitensor.a $(OBJECTS)
	@ranlib $(ITENSOR_LIBDIR)/libitensor.a

libitensor-g.a: mkdebugdir $(GOBJECTS)
	@echo "Building static library $(ITENSOR_LIBDIR)/libitensor-g.a"
	@ar r $(ITENSOR_LIBDIR)/libitensor-g.a $(GOBJECTS)
	@ranlib $(ITENSOR_LIBDIR)/libitensor-g.a

build: libitensor.a

debug: libitensor-g.a

mkdebugdir:
	@mkdir -p .debug_objs
	@mkdir -p .debug_objs/util
	@mkdir -p .debug_objs/tensor
	@mkdir -p .debug_objs/itdata
	@mkdir -p .debug_objs/mps

clean:	
	@rm -fr *.o .debug_objs util/*.o tensor/*.o itdata/*.o mps/*.o libitensor.a libitensor-g.a

GDEPHEADERS=real.h global.h index.h util/readwrite.h
GDEPHEADERS+= tensor/types.h tensor/vecrange.h tensor/ten.h tensor/ten.ih \
tensor/teniter.h tensor/range.h tensor/lapack_wrap.h tensor/vec.h
tensor/vec.o: $(GDEPHEADERS)
.debug_objs/tensor/vec.o: $(GDEPHEADERS)
GDEPHEADERS+= tensor/matrange.h  tensor/mat.h
tensor/mat.o: $(GDEPHEADERS)
.debug_objs/tensor/mat.o: $(GDEPHEADERS)
GDEPHEADERS+= tensor/slicemat.h tensor/algs.h
tensor/algs.o: $(GDEPHEADERS)
.debug_objs/tensor/algs.o: $(GDEPHEADERS)
GDEPHEADERS+= tensor/permutation.h tensor/slicerange.h tensor/sliceten.h \
tensor/contract.h itdata/task_types.h indexset.ih indexset.h
tensor/contract.o: $(GDEPHEADERS)
.debug_objs/tensor/contract.o: $(GDEPHEADERS)
ITDEPHEADERS= itdata/dense.h
itdata/dense.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itdata/dense.o: $(ITDEPHEADERS) $(GDEPHEADERS)
ITDEPHEADERS+= itdata/diag.h
itdata/diag.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itdata/diag.o: $(ITDEPHEADERS) $(GDEPHEADERS)
ITDEPHEADERS+= itdata/itcombiner.h
itdata/itcombiner.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itdata/itcombiner.o: $(ITDEPHEADERS) $(GDEPHEADERS)
ITDEPHEADERS+= itdata/iqtreal.h
itdata/iqtreal.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itdata/iqtreal.o: $(ITDEPHEADERS) $(GDEPHEADERS)
ITDEPHEADERS+= itdata/iqtcombiner.h
itdata/iqtcombiner.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itdata/iqtcombiner.o: $(ITDEPHEADERS) $(GDEPHEADERS)
ITDEPHEADERS+= itdata/iqtdiag.h
itdata/iqtdiag.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itdata/iqtdiag.o: $(ITDEPHEADERS) $(GDEPHEADERS)
ITDEPHEADERS+= index.h
index.o: $(ITDEPHEADERS)
.debug_objs/index.o: $(ITDEPHEADERS)
ITDEPHEADERS+= indexset.ih indexset.h
indexset.o: $(ITDEPHEADERS)
.debug_objs/indexset.o: $(ITDEPHEADERS)
ITDEPHEADERS+= itensor_interface.h itensor_interface.ih \
itensor.ih itensor.h itdata/itdata.h itdata/dense.h itdata/diag.h
itensor.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itensor.o: $(ITDEPHEADERS) $(GDEPHEADERS)
ITDEPHEADERS+= qn.h iqindex.h
iqindex.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/iqindex.o: $(ITDEPHEADERS) $(GDEPHEADERS)
ITDEPHEADERS+= iqtensor.ih iqtensor.h detail/skip_iterator.h
iqtensor.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/iqtensor.o: $(ITDEPHEADERS) $(GDEPHEADERS)
GDEPHEADERS+= spectrum.h
spectrum.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/spectrum.o: $(ITDEPHEADERS) $(GDEPHEADERS)
GDEPHEADERS+= svdalgs.h
svdalgs.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/svdalgs.o: $(ITDEPHEADERS) $(GDEPHEADERS)
GDEPHEADERS+= mps/mps.h
mps/mps.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/mps/mps.o: $(ITDEPHEADERS) $(GDEPHEADERS)
GDEPHEADERS+= mps/mpo.h
mps/mpo.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/mps/mpo.o: $(ITDEPHEADERS) $(GDEPHEADERS)
