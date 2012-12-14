include ../../this_dir.mk
include ../../options.mk
################################################################

LIBNAMES=itensor matrix utilities

#Define Dirs -----------
INCLUDEDIR=$(PREFIX)/include
INCLUDEFLAGS=. -I$(INCLUDEDIR) -I$(BOOST_DIR) $(BLAS_LAPACK_INCLUDEFLAGS)

#Mappings --------------
LOCAL_LIBFLAGS=$(patsubst %,-l%, $(LIBNAMES))
LOCAL_LIBGFLAGS=$(patsubst %,-l%-g, $(LIBNAMES))
LIBFILES=$(patsubst %,$(LIBDIR)/lib%.a, $(LIBNAMES))
LIBGFILES=$(patsubst %,$(LIBDIR)/lib%-g.a, $(LIBNAMES))
REL_TENSOR_HEADERS=$(patsubst %,$(INCLUDEDIR)/%, $(TENSOR_HEADERS))

#Define Flags ----------
CCFLAGS=$(CPPFLAGS) -I$(INCLUDEFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= -I$(INCLUDEFLAGS) -DDEBUG -DMATRIXBOUNDS -DITENSOR_USE_AT -DBOUNDS -g -Wall -ansi
LIBFLAGS=-L$(LIBDIR) $(LOCAL_LIBFLAGS) $(BLAS_LAPACK_LIBFLAGS)
LIBGFLAGS=-L$(LIBDIR) $(LOCAL_LIBGFLAGS) $(BLAS_LAPACK_LIBFLAGS)

#Rules ------------------

%.o: %.cc $(LIBFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(LIBFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: one

debug: one-g

one: one.o $(LIBFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) one.o -o one $(LIBFLAGS)

one-g: mkdebugdir .debug_objs/one.o $(LIBGFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/one.o -o one-g $(LIBGFLAGS)

mkdebugdir:
	mkdir -p .debug_objs

clean:
	rm -fr *.o .debug_objs one one-g
