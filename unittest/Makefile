include ../this_dir.mk
include ../options.mk
include unittest_options.mk
################################################################

REPORT_LEVEL=confirm
#REPORT_LEVEL=detailed

HEADERS=test.h
SOURCES=test.cc matrix_test.cc real_test.cc index_test.cc itensor_test.cc\
combiner_test.cc iqtensor_test.cc mps_test.cc

TENSOR_HEADERS=dmrg.h hams.h
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
REL_TENSOR_HEADERS=$(patsubst %,$(INCLUDEDIR)/%, $(TENSOR_HEADERS))
OBJECTS=$(patsubst %.cc,%.o, $(SOURCES))
GOBJECTS=$(patsubst %.cc,.debug_objs/%.o, $(SOURCES))

#Define Flags ----------
CCFLAGS=$(CPPFLAGS) -I$(INCLUDEFLAGS) -O2 
CCGFLAGS= -I$(INCLUDEFLAGS) -DDEBUG -DMATRIXBOUNDS -DBOUNDS -g -Wall -ansi
LIBFLAGS=-L$(THIS_LIBDIR) $(LOCAL_LIBFLAGS) $(BLAS_LAPACK_LIBFLAGS) $(BOOST_UNITTEST_LIBFLAGS)
LIBGFLAGS=-L$(THIS_LIBDIR) $(LOCAL_LIBGFLAGS) $(BLAS_LAPACK_LIBFLAGS) $(BOOST_UNITTEST_LIBFLAGS)

#Rules ------------------

%.o: %.cc $(LIBFILES) $(REL_TENSOR_HEADERS) $(HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(LIBFILES) $(REL_TENSOR_HEADERS) $(HEADERS)
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
