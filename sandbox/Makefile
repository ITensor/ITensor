include ../this_dir.mk
include ../options.mk
################################################################

TENSOR_HEADERS=dmrg.h hams.h
LIBNAMES=matrix utilities itensor

#################################################################

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
CCFLAGS=$(CPPFLAGS) -I$(INCLUDEFLAGS) -O2
CCGFLAGS= -I$(INCLUDEFLAGS) -DDEBUG -DMATRIXBOUNDS -DITENSOR_USE_AT -DBOUNDS -g -Wall -ansi
LIBFLAGS=-L$(LIBDIR) $(LOCAL_LIBFLAGS) $(BLAS_LAPACK_LIBFLAGS)
LIBGFLAGS=-L$(LIBDIR) $(LOCAL_LIBGFLAGS) $(BLAS_LAPACK_LIBFLAGS)

#Rules ------------------

%.o: %.cc $(LIBFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(LIBFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: dmrg iqdmrg onesiteopt iqonesiteopt

debug: dmrg-g iqdmrg-g onesiteopt-g iqonesiteopt

dmrg: dmrg.o $(LIBFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) dmrg.o -o dmrg $(LIBFLAGS)

dmrg-g: mkdebugdir .debug_objs/dmrg.o $(LIBGFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/dmrg.o -o dmrg-g $(LIBGFLAGS)

iqdmrg: iqdmrg.o $(LIBFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) iqdmrg.o -o iqdmrg $(LIBFLAGS)

iqdmrg-g: mkdebugdir .debug_objs/iqdmrg.o $(LIBGFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/iqdmrg.o -o iqdmrg-g $(LIBGFLAGS)

onesiteopt: onesiteopt.o $(LIBFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) onesiteopt.o -o onesiteopt $(LIBFLAGS)

onesiteopt-g: mkdebugdir .debug_objs/onesiteopt.o $(LIBGFILES) $(REL_TENSOR_HEADERS) 
	$(CCCOM) $(CCGFLAGS) .debug_objs/onesiteopt.o -o onesiteopt-g $(LIBGFLAGS)

iqonesiteopt: iqonesiteopt.o $(LIBFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) iqonesiteopt.o -o iqonesiteopt $(LIBFLAGS)

iqonesiteopt-g: mkdebugdir .debug_objs/iqonesiteopt.o $(LIBGFILES) $(REL_TENSOR_HEADERS) 
	$(CCCOM) $(CCGFLAGS) .debug_objs/iqonesiteopt.o -o iqonesiteopt-g $(LIBGFLAGS)


mkdebugdir:
	mkdir -p .debug_objs

clean:
	rm -fr *.o .debug_objs dmrg dmrg-g iqdmrg iqdmrg-g onesiteopt onesiteopt-g iqonesiteopt iqonesiteopt-g
