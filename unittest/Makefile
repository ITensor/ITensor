include ../this_dir.mk
include ../options.mk
################################################################

HEADERS=test.h
SOURCES = test.cc
#SOURCES+= matrix_test.cc
#SOURCES+= real_test.cc
#SOURCES+= option_test.cc
#SOURCES+= index_test.cc
#SOURCES+= indexset_test.cc
#SOURCES+= tensor_test.cc
#SOURCES+= contract_test.cc
SOURCES+= itensor_test.cc
#SOURCES+= combiner_test.cc
SOURCES+= iqindex_test.cc
#SOURCES+= iqcombiner_test.cc
#SOURCES+= iqtensor_test.cc
#SOURCES+= mps_test.cc
#SOURCES+= mpo_test.cc
#SOURCES+= eigensolver_test.cc
#SOURCES+= regression_test.cc
#SOURCES+= spectrum_test.cc
#SOURCES+= svdalgs_test.cc
#SOURCES+= webpage_test.cc
#SOURCES+= localmpo_test.cc
#SOURCES+= siteset_test.cc
#SOURCES+= bondgate_test.cc

##################################################################

OBJECTS=$(patsubst %.cc,%.o, $(SOURCES))
GOBJECTS=$(patsubst %.cc,.debug_objs/%.o, $(SOURCES))

#Define Flags ----------
CCFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(CPPFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS)
LIBFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBFLAGS)
LIBGFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBGFLAGS)

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
	@./test-g

growl:
	@make run && growlnotify -m "Tests Passed" || growlnotify -m "Test Failed"

debug: test-g

#build: test


#test: test.o $(ITENSOR_LIBS) $(OBJECTS)
#	$(CCCOM) $(CCFLAGS) $(OBJECTS) -o test $(LIBFLAGS)

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

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/combiner.h
combiner_test.o: $(LIBHEADERS)
.debug_objs/combiner_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/iqcombiner.h
iqcombiner_test.o: $(LIBHEADERS)
.debug_objs/iqcombiner_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/mps.h
mps_test.o: $(LIBHEADERS)
.debug_objs/mps_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/mpo.h
mpo_test.o: $(LIBHEADERS)
.debug_objs/mpo_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/svdalgs.h
svdalgs_test.o: $(LIBHEADERS)
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

