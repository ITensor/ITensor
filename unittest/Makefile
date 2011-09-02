include ../this_dir.mk
include ../options.mk
include unittest_options.mk
################################################################

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

#Define Flags ----------
CCFLAGS=$(CPPFLAGS) -I$(INCLUDEFLAGS) -O2  -DNDEBUG -DBOOST_DISABLE_ASSERTS
CCGFLAGS= -I$(INCLUDEFLAGS) -DMATRIXBOUNDS -DITENSOR_USE_AT -DBOUNDS -g -Wall -ansi
LIBFLAGS=-L$(THIS_LIBDIR) $(LOCAL_LIBFLAGS) $(BLAS_LAPACK_LIBFLAGS) $(BOOST_UNITTEST_LIBFLAGS)
LIBGFLAGS=-L$(THIS_LIBDIR) $(LOCAL_LIBGFLAGS) $(BLAS_LAPACK_LIBFLAGS) $(BOOST_UNITTEST_LIBFLAGS)

#Rules ------------------

%.o: %.cc $(LIBFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(LIBFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: hello

debug: hello-g

hello: hello.o $(LIBFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) hello.o -o hello $(LIBFLAGS)

hello-g: mkdebugdir .debug_objs/hello.o $(LIBGFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/hello.o -o hello-g $(LIBGFLAGS)


mkdebugdir:
	mkdir -p .debug_objs

clean:
	rm -fr *.o .debug_objs 
