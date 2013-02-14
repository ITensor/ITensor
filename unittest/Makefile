include ../this_dir.mk
include ../options.mk
include unittest_options.mk
################################################################

REPORT_LEVEL=confirm
#REPORT_LEVEL=detailed

HEADERS=test.h
SOURCES = test.cc
#SOURCES+= matrix_test.cc
SOURCES+= real_test.cc
SOURCES+= index_test.cc
SOURCES+= itensor_test.cc
#SOURCES+= itsparse_test.cc
#SOURCES+= combiner_test.cc
SOURCES+= iqcombiner_test.cc
SOURCES+= iqtensor_test.cc
#SOURCES+= mps_test.cc
#SOURCES+= mpo_test.cc
#SOURCES+= eigensolver_test.cc
#SOURCES+= regression_test.cc
SOURCES+= svdworker_test.cc
#SOURCES+= iqtsparse_test.cc
#SOURCES+= webpage_test.cc
#SOURCES+= localmpo_test.cc
SOURCES+= option_test.cc
SOURCES+= indexset_test.cc

LIBNAMES=matrix utilities itensor

#################################################################

#Define Dirs -----------
INCLUDEDIR=$(PREFIX)/include
INCLUDEFLAGS=. -I$(INCLUDEDIR) -I$(BOOST_DIR) $(BLAS_LAPACK_INCLUDEFLAGS) -I$(BOOST_UNITTEST_INCLUDEDIR)
THIS_LIBDIR=$(LIBDIR) -L$(BOOST_UNITTEST_LIBDIR)

#Mappings --------------
LOCAL_LIBFLAGS=$(patsubst %,-l%, $(LIBNAMES))
LOCAL_LIBGFLAGS=$(patsubst %,-l%-g, $(LIBNAMES))
LIBFILES=$(patsubst %,$(LIBDIR)/lib%.a, $(LIBNAMES))
LIBGFILES=$(patsubst %,$(LIBDIR)/lib%-g.a, $(LIBNAMES))
OBJECTS=$(patsubst %.cc,%.o, $(SOURCES))
GOBJECTS=$(patsubst %.cc,.debug_objs/%.o, $(SOURCES))

#Define Flags ----------
CCFLAGS=$(CPPFLAGS) -I$(INCLUDEFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= -I$(INCLUDEFLAGS) $(DEBUGFLAGS)
LIBFLAGS=-L$(THIS_LIBDIR) $(LOCAL_LIBFLAGS) $(BLAS_LAPACK_LIBFLAGS) $(BOOST_UNITTEST_LIBFLAGS)
LIBGFLAGS=-L$(THIS_LIBDIR) $(LOCAL_LIBGFLAGS) $(BLAS_LAPACK_LIBFLAGS) $(BOOST_UNITTEST_LIBFLAGS)

#Rules ------------------

%.o: %.cc $(HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

run: test-g
	@echo 
	@echo Running all tests...
	@echo 
	@./test-g --report_level=$(REPORT_LEVEL)

growl:
	@make run && growlnotify -m "Tests Passed" || growlnotify -m "Test Failed"

build: test

debug: test-g

test: test.o $(LIBFILES) $(REL_TENSOR_HEADERS) $(HEADERS) $(OBJECTS)
	$(CCCOM) $(CCFLAGS) $(OBJECTS) -o test $(LIBFLAGS)

test-g: mkdebugdir .debug_objs/test.o $(LIBGFILES) $(REL_TENSOR_HEADERS) $(HEADERS) $(GOBJECTS)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o test-g $(LIBGFLAGS)


mkdebugdir:
	mkdir -p .debug_objs

clean:
	rm -fr *.o .debug_objs test test-g

LIBHEADERS=$(INCLUDEDIR)/matrix.h
matrix_test.o: $(LIBHEADERS)
.debug/matrix_test.o: $(LIBHEADERS)

LIBHEADERS=$(INCLUDEDIR)/real.h
real_test.o: $(LIBHEADERS)
.debug/real_test.o: $(LIBHEADERS)

LIBHEADERS+= $(INCLUDEDIR)/index.h
index_test.o: $(LIBHEADERS)
.debug_objs/index_test.o: $(LIBHEADERS)

LIBHEADERS+= $(INCLUDEDIR)/itensor.h
itensor_test.o: $(LIBHEADERS)
.debug_objs/itensor_test.o: $(LIBHEADERS)

LIBHEADERS+= $(INCLUDEDIR)/itsparse.h
itsparse_test.o: $(LIBHEADERS)
.debug_objs/itsparse_test.o: $(LIBHEADERS)

LIBHEADERS+= $(INCLUDEDIR)/combiner.h
combiner_test.o: $(LIBHEADERS)
.debug_objs/combiner_test.o: $(LIBHEADERS)

LIBHEADERS+= $(INCLUDEDIR)/iqcombiner.h
iqcombiner_test.o: $(LIBHEADERS)
.debug_objs/iqcombiner_test.o: $(LIBHEADERS)

LIBHEADERS+= $(INCLUDEDIR)/iqtsparse.h
iqtensor_test.o: $(LIBHEADERS)
.debug_objs/iqtensor_test.o: $(LIBHEADERS)

LIBHEADERS+= $(INCLUDEDIR)/iqtsparse.h
iqtsparse_test.o: $(LIBHEADERS)
.debug_objs/iqtsparse_test.o: $(LIBHEADERS)

LIBHEADERS+= $(INCLUDEDIR)/mps.h
mps_test.o: $(LIBHEADERS)
.debug_objs/mps_test.o: $(LIBHEADERS)

LIBHEADERS+= $(INCLUDEDIR)/mpo.h
mpo_test.o: $(LIBHEADERS)
.debug_objs/mpo_test.o: $(LIBHEADERS)

LIBHEADERS+= $(INCLUDEDIR)/svdworker.h
svdworker_test.o: $(LIBHEADERS)
.debug_objs/svdworker_test.o: $(LIBHEADERS)

regression_test.o: $(LIBHEADERS)
.debug_objs/regression_test.o: $(LIBHEADERS)

eigensolver_test.o: $(INCLUDEDIR)/eigensolver.h
.debug_objs/eigensolver_test.o: $(INCLUDEDIR)/eigensolver.h

webpage_test.o: $(INCLUDEDIR)/iqtensor.h
.debug_objs/webpage_test.o: $(INCLUDEDIR)/iqtensor.h

LIBHEADERS+= $(INCLUDEDIR)/localmpo.h
localmpo_test.o: $(LIBHEADERS)
.debug_objs/localmpo_test.o: $(LIBHEADERS)

