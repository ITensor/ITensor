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
SOURCES+= itdata/combiner.cc 
SOURCES+= itdata/diag.cc 
SOURCES+= itdata/qdense.cc
SOURCES+= itdata/qcombiner.cc 
SOURCES+= itdata/qdiag.cc
SOURCES+= itdata/qmixed.cc
SOURCES+= itdata/scalar.cc 
##SOURCES+= itdata/itlazy.cc
SOURCES+= index.cc 
SOURCES+= itensor_interface.cc 
SOURCES+= itensor_operators.cc 
SOURCES+= itensor.cc 
SOURCES+= qn.cc 
SOURCES+= iqindex.cc 
SOURCES+= iqtensor.cc 
SOURCES+= spectrum.cc 
SOURCES+= decomp.cc 
SOURCES+= svd.cc 
SOURCES+= hermitian.cc 
SOURCES+= global.cc
SOURCES+= mps/mps.cc 
SOURCES+= mps/mpsalgs.cc 
SOURCES+= mps/mpo.cc 
SOURCES+= mps/mpoalgs.cc 
SOURCES+= mps/autompo.cc

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

$(ITENSOR_LIBDIR)/libitensor.a: $(OBJECTS)
	@echo "Building static library $(ITENSOR_LIBDIR)/libitensor.a"
	@ar r $(ITENSOR_LIBDIR)/libitensor.a $(OBJECTS)
	@ranlib $(ITENSOR_LIBDIR)/libitensor.a

$(ITENSOR_LIBDIR)/libitensor-g.a: mkdebugdir $(GOBJECTS)
	@echo "Building static library $(ITENSOR_LIBDIR)/libitensor-g.a"
	@ar r $(ITENSOR_LIBDIR)/libitensor-g.a $(GOBJECTS)
	@ranlib $(ITENSOR_LIBDIR)/libitensor-g.a

$(ITENSOR_LIBDIR)/libitensor.$(DYLIB_EXT): $(OBJECTS)
ifdef ITENSOR_MAKE_DYLIB
ifeq ($(ITENSOR_MAKE_DYLIB),1)
	@echo "Building dynamic library $(ITENSOR_LIBDIR)/libitensor.$(DYLIB_EXT)"
	@$(CCCOM) $(DYLIB_FLAGS) -o $(ITENSOR_LIBDIR)/libitensor.$(DYLIB_EXT) $(OBJECTS) $(BLAS_LAPACK_LIBFLAGS)
endif
endif

$(ITENSOR_LIBDIR)/libitensor-g.$(DYLIB_EXT): mkdebugdir $(GOBJECTS)
ifdef ITENSOR_MAKE_DYLIB
ifeq ($(ITENSOR_MAKE_DYLIB),1)
	@echo "Building dynamic library $(ITENSOR_LIBDIR)/libitensor-g.$(DYLIB_EXT)"
	@$(CCCOM) $(DYLIB_FLAGS) -o $(ITENSOR_LIBDIR)/libitensor-g.$(DYLIB_EXT) $(GOBJECTS) $(BLAS_LAPACK_LIBFLAGS)
endif
endif

touch_all_headers:
	@touch all.h
	@touch all_basic.h
	@touch all_mps.h

build: $(ITENSOR_LIBDIR)/libitensor.a $(ITENSOR_LIBDIR)/libitensor.$(DYLIB_EXT) touch_all_headers

debug: $(ITENSOR_LIBDIR)/libitensor-g.a $(ITENSOR_LIBDIR)/libitensor-g.$(DYLIB_EXT) touch_all_headers

mkdebugdir:
	@mkdir -p .debug_objs
	@mkdir -p .debug_objs/util
	@mkdir -p .debug_objs/tensor
	@mkdir -p .debug_objs/itdata
	@mkdir -p .debug_objs/mps

clean:	
	@rm -fr *.o .debug_objs util/*.o tensor/*.o itdata/*.o mps/*.o libitensor.a libitensor-g.a

util/input.o: util/input.h
.debug_objs/util/input.o: util/input.h

GDEPHEADERS=real.h global.h index.h util/readwrite.h
GDEPHEADERS+= tensor/types.h tensor/vecrange.h tensor/ten.h tensor/ten_impl.h \
tensor/teniter.h tensor/range.h tensor/lapack_wrap.h tensor/vec.h util/safe_ptr.h
tensor/vec.o: $(GDEPHEADERS)
.debug_objs/tensor/vec.o: $(GDEPHEADERS)
GDEPHEADERS+= tensor/matrange.h  tensor/mat.h
tensor/mat.o: $(GDEPHEADERS)
.debug_objs/tensor/mat.o: $(GDEPHEADERS)
tensor/gemm.o: $(GDEPHEADERS)
.debug_objs/tensor/gemm.o: $(GDEPHEADERS)
GDEPHEADERS+= tensor/slicemat.h tensor/algs.h tensor/algs_impl.h
tensor/algs.o: $(GDEPHEADERS)
.debug_objs/tensor/algs.o: $(GDEPHEADERS)
GDEPHEADERS+= tensor/permutation.h tensor/slicerange.h tensor/sliceten.h \
tensor/contract.h itdata/task_types.h indexset_impl.h indexset.h
tensor/contract.o: $(GDEPHEADERS)
.debug_objs/tensor/contract.o: $(GDEPHEADERS)
ITDEPHEADERS= itdata/dense.h 
itdata/dense.o: $(ITDEPHEADERS) $(GDEPHEADERS) util/tensorstats.h
.debug_objs/itdata/dense.o: $(ITDEPHEADERS) $(GDEPHEADERS) util/tensorstats.h
ITDEPHEADERS+= itdata/diag.h
itdata/diag.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itdata/diag.o: $(ITDEPHEADERS) $(GDEPHEADERS)
ITDEPHEADERS+= itdata/combiner.h
itdata/combiner.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itdata/combiner.o: $(ITDEPHEADERS) $(GDEPHEADERS)
ITDEPHEADERS+= itdata/qdense.h itdata/qutil.h
itdata/qdense.o: $(ITDEPHEADERS) $(GDEPHEADERS) util/tensorstats.h
.debug_objs/itdata/qdense.o: $(ITDEPHEADERS) $(GDEPHEADERS) util/tensorstats.h
ITDEPHEADERS+= itdata/qcombiner.h
itdata/qcombiner.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itdata/qcombiner.o: $(ITDEPHEADERS) $(GDEPHEADERS)
ITDEPHEADERS+= itdata/qdiag.h
itdata/qdiag.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itdata/qdiag.o: $(ITDEPHEADERS) $(GDEPHEADERS)
ITDEPHEADERS+= itdata/qmixed.h
itdata/qmixed.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itdata/qmixed.o: $(ITDEPHEADERS) $(GDEPHEADERS)
ITDEPHEADERS+= itdata/scalar.h
itdata/scalar.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itdata/scalar.o: $(ITDEPHEADERS) $(GDEPHEADERS)
ITDEPHEADERS+= index.h
index.o: $(ITDEPHEADERS)
.debug_objs/index.o: $(ITDEPHEADERS)
ITDEPHEADERS+= indexset_impl.h indexset.h
indexset.o: $(ITDEPHEADERS)
.debug_objs/indexset.o: $(ITDEPHEADERS)
ITDEPHEADERS+= itensor_interface.h itensor_interface_impl.h \
itensor_impl.h itensor.h itdata/itdata.h itdata/dense.h itdata/diag.h
itensor_interface.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itensor_interface.o: $(ITDEPHEADERS) $(GDEPHEADERS)
itensor_operators.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itensor_operators.o: $(ITDEPHEADERS) $(GDEPHEADERS)
itensor.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itensor.o: $(ITDEPHEADERS) $(GDEPHEADERS)
ITDEPHEADERS+= qn.h
qn.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/qn.o: $(ITDEPHEADERS) $(GDEPHEADERS)
ITDEPHEADERS+= iqindex.h
iqindex.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/iqindex.o: $(ITDEPHEADERS) $(GDEPHEADERS)
ITDEPHEADERS+= iqtensor_impl.h iqtensor.h detail/skip_iterator.h
iqtensor.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/iqtensor.o: $(ITDEPHEADERS) $(GDEPHEADERS)
GDEPHEADERS+= spectrum.h
spectrum.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/spectrum.o: $(ITDEPHEADERS) $(GDEPHEADERS)
GDEPHEADERS+= decomp.h
decomp.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/decomp.o: $(ITDEPHEADERS) $(GDEPHEADERS)
svd.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/svd.o: $(ITDEPHEADERS) $(GDEPHEADERS)
hermitian.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/hermitian.o: $(ITDEPHEADERS) $(GDEPHEADERS)
GDEPHEADERS+= mps/mps.h
mps/mps.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/mps/mps.o: $(ITDEPHEADERS) $(GDEPHEADERS)
mps/mpsalgs.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/mps/mpsalgs.o: $(ITDEPHEADERS) $(GDEPHEADERS)
GDEPHEADERS+= mps/mpo.h
mps/mpo.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/mps/mpo.o: $(ITDEPHEADERS) $(GDEPHEADERS)
mps/mpoalgs.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/mps/mpoalgs.o: $(ITDEPHEADERS) $(GDEPHEADERS)
mps/autompo.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/mps/autompo.o: $(ITDEPHEADERS) $(GDEPHEADERS)
