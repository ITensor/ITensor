include ../this_dir.mk
include ../options.mk
include unittest_options.mk
################################################################

REPORT_LEVEL=confirm
#REPORT_LEVEL=detailed

HEADERS=test.h
SOURCES = test.cc
SOURCES+= matrix_test.cc
SOURCES+= real_test.cc
SOURCES+= index_test.cc
SOURCES+= itensor_test.cc
SOURCES+= itsparse_test.cc
SOURCES+= combiner_test.cc
SOURCES+= iqindex_test.cc
SOURCES+= iqcombiner_test.cc
SOURCES+= iqtensor_test.cc
SOURCES+= mps_test.cc
SOURCES+= mpo_test.cc
SOURCES+= eigensolver_test.cc
SOURCES+= regression_test.cc
SOURCES+= svdworker_test.cc
SOURCES+= iqtsparse_test.cc
SOURCES+= webpage_test.cc
SOURCES+= localmpo_test.cc
SOURCES+= option_test.cc
SOURCES+= indexset_test.cc

##################################################################

OBJECTS=$(patsubst %.cc,%.o, $(SOURCES))
GOBJECTS=$(patsubst %.cc,.debug_objs/%.o, $(SOURCES))

#Define Flags ----------
CCFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(CPPFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) -I$(BOOST_UNITTEST_INCLUDEDIR) $(DEBUGFLAGS)
LIBFLAGS=-L$(ITENSOR_LIBDIR) -L$(BOOST_UNITTEST_LIBDIR) $(BOOST_UNITTEST_LIBFLAGS) $(ITENSOR_LIBFLAGS)
LIBGFLAGS=-L$(ITENSOR_LIBDIR) -L$(BOOST_UNITTEST_LIBDIR) $(BOOST_UNITTEST_LIBFLAGS) $(ITENSOR_LIBGFLAGS)

#Rules ------------------

%.o: %.cc $(ITENSOR_LIBS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(ITENSOR_GLIBS)
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

test: test.o $(ITENSOR_LIBS) $(OBJECTS)
	$(CCCOM) $(CCFLAGS) $(OBJECTS) -o test $(LIBFLAGS)

test-g: mkdebugdir .debug_objs/test.o $(ITENSOR_GLIBS) $(GOBJECTS)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o test-g $(LIBGFLAGS)


mkdebugdir:
	mkdir -p .debug_objs

clean:
	rm -fr *.o .debug_objs test test-g

LIBHEADERS=$(ITENSOR_INCLUDEDIR)/matrix.h
matrix_test.o: $(LIBHEADERS)
.debug/matrix_test.o: $(LIBHEADERS)

LIBHEADERS=$(ITENSOR_INCLUDEDIR)/real.h
real_test.o: $(LIBHEADERS)
.debug/real_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/index.h
index_test.o: $(LIBHEADERS)
.debug_objs/index_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/itensor.h
itensor_test.o: $(LIBHEADERS)
.debug_objs/itensor_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/itsparse.h
itsparse_test.o: $(LIBHEADERS)
.debug_objs/itsparse_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/combiner.h
combiner_test.o: $(LIBHEADERS)
.debug_objs/combiner_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/iqcombiner.h
iqcombiner_test.o: $(LIBHEADERS)
.debug_objs/iqcombiner_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/iqtsparse.h
iqtensor_test.o: $(LIBHEADERS)
.debug_objs/iqtensor_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/iqtsparse.h
iqtsparse_test.o: $(LIBHEADERS)
.debug_objs/iqtsparse_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/mps.h
mps_test.o: $(LIBHEADERS)
.debug_objs/mps_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/mpo.h
mpo_test.o: $(LIBHEADERS)
.debug_objs/mpo_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/svdworker.h
svdworker_test.o: $(LIBHEADERS)
.debug_objs/svdworker_test.o: $(LIBHEADERS)

regression_test.o: $(LIBHEADERS)
.debug_objs/regression_test.o: $(LIBHEADERS)

eigensolver_test.o: $(ITENSOR_INCLUDEDIR)/eigensolver.h
.debug_objs/eigensolver_test.o: $(ITENSOR_INCLUDEDIR)/eigensolver.h

webpage_test.o: $(ITENSOR_INCLUDEDIR)/iqtensor.h
.debug_objs/webpage_test.o: $(ITENSOR_INCLUDEDIR)/iqtensor.h

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/localmpo.h
localmpo_test.o: $(LIBHEADERS)
.debug_objs/localmpo_test.o: $(LIBHEADERS)

