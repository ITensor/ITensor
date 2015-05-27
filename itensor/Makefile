#
# Makefile for ITensor library
#

include ../this_dir.mk
include ../options.mk

####################################

SOURCES+= util/args.cc     
SOURCES+= util/input.cc
SOURCES+= matrix/vec.cc 
SOURCES+= matrix/mat.cc 
SOURCES+= matrix/algs.cc 
SOURCES+= tensor/contract.cc 
SOURCES+= index.cc 
SOURCES+= itensor.cc 
SOURCES+= iqindex.cc 
SOURCES+= itdata/iqtdata.cc
SOURCES+= iqtensor.cc 
SOURCES+= spectrum.cc 
SOURCES+= svdalgs.cc 
SOURCES+= mps/mps.cc 
#SOURCES+= mps/mpo.cc 
#SOURCES+= mps/tevol.cc
#SOURCES+= mps/autompo.cc

####################################

CCFLAGS= $(ITENSOR_INCLUDEFLAGS) $(OPTIMIZATIONS) -DPLATFORM_$(PLATFORM) -D__ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES=0
CCGFLAGS= $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS) -DPLATFORM_$(PLATFORM) -D__ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES=0

OBJECTS= $(patsubst %.cc,%.o, $(SOURCES))
GOBJECTS= $(patsubst %.cc,.debug_objs/%.o, $(SOURCES))

.SUFFIXES: .cc $(SUFFIXES)

%.o: %.cc
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

install: installo installg
	@echo $(OBJECTS)

libitensor.a: $(OBJECTS)
	ar rv libitensor.a $(OBJECTS)
	ranlib libitensor.a

libitensor-g.a: mkdebugdir $(GOBJECTS)
	ar rv libitensor-g.a $(GOBJECTS)
	ranlib libitensor-g.a

libs: libitensor.a 
	cp libitensor.a $(ITENSOR_LIBDIR)

installo: libitensor.a
	cp libitensor.a $(ITENSOR_LIBDIR)
	ranlib $(ITENSOR_LIBDIR)/libitensor.a

debug: installg
installg: libitensor-g.a
	cp libitensor-g.a $(ITENSOR_LIBDIR)
	ranlib $(ITENSOR_LIBDIR)/libitensor-g.a

mkdebugdir:
	mkdir -p .debug_objs
	mkdir -p .debug_objs/util
	mkdir -p .debug_objs/matrix
	mkdir -p .debug_objs/tensor
	mkdir -p .debug_objs/itdata
	mkdir -p .debug_objs/mps

clean:	
	rm -fr *.o .debug_objs util/*.o matrix/*.o tensor/*.o itdata/*.o mps/*.o libitensor.a libitensor-g.a

GDEPHEADERS=real.h global.h index.h
GDEPHEADERS+= matrix/strideiter.h matrix/vec.h
matrix/vec.o: $(GDEPHEADERS)
.debug_objs/matrix/vec.o: $(GDEPHEADERS)
GDEPHEADERS+= matrix/mrange.h matrix/matiter.h matrix/mat.h
matrix/mat.o: $(GDEPHEADERS)
.debug_objs/matrix/mat.o: $(GDEPHEADERS)
GDEPHEADERS+= matrix/slicemat.h matrix/algs.h
matrix/algs.o: $(GDEPHEADERS)
.debug_objs/matrix/algs.o: $(GDEPHEADERS)
GDEPHEADERS+= tensor/range.h tensor/ten.h tensor/permutation.h \
tensor/permute.h tensor/contract.h
tensor/contract.o: $(GDEPHEADERS)
.debug_objs/tensor/contract.o: $(GDEPHEADERS)
ITDEPHEADERS= index.h
index.o: $(ITDEPHEADERS)
.debug_objs/index.o: $(ITDEPHEADERS)
ITDEPHEADERS+= indexset.ih indexset.h
indexset.o: $(ITDEPHEADERS)
.debug_objs/indexset.o: $(ITDEPHEADERS)
ITDEPHEADERS+= task_types.h itensor.ih itensor.h itdata/itdata.h itdata/itreal.h itdata/itdiag.h itdata/itcplx.h
itensor.o: $(ITDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itensor.o: $(ITDEPHEADERS) $(GDEPHEADERS)
IQDEPHEADERS= qn.h iqindex.h
iqindex.o: $(IQDEPHEADERS) $(GDEPHEADERS)
.debug_objs/iqindex.o: $(IQDEPHEADERS) $(GDEPHEADERS)
IQDEPHEADERS+= $(ITDEPHEADERS) itdata/iqtdata.h
itdata/iqtdata.o: $(IQDEPHEADERS) $(GDEPHEADERS)
.debug_objs/itdata/iqtdata.o: $(IQDEPHEADERS) $(GDEPHEADERS)
IQDEPHEADERS+= iqtensor.ih iqtensor.h detail/skip_iterator.h
iqtensor.o: $(IQDEPHEADERS) $(GDEPHEADERS)
.debug_objs/iqtensor.o: $(IQDEPHEADERS) $(GDEPHEADERS)
GDEPHEADERS+= spectrum.h
spectrum.o: $(IQDEPHEADERS) $(GDEPHEADERS)
.debug_objs/spectrum.o: $(IQDEPHEADERS) $(GDEPHEADERS)
GDEPHEADERS+= svdalgs.h
svdalgs.o: $(IQDEPHEADERS) $(GDEPHEADERS)
.debug_objs/svdalgs.o: $(IQDEPHEADERS) $(GDEPHEADERS)
