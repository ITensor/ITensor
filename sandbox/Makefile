include ../this_dir.mk
include ../options.mk
################################################################

TENSOR_HEADERS=core.h

#################################################################

#Mappings --------------
REL_TENSOR_HEADERS=$(patsubst %,$(ITENSOR_INCLUDEDIR)/%, $(TENSOR_HEADERS))

#Define Flags ----------
CCFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(CPPFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS)
LIBFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBFLAGS)
LIBGFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBGFLAGS)

#Rules ------------------

%.o: %.cc $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(ITENSOR_GLIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: dmrg iqdmrg

debug: dmrg-g iqdmrg-g

all: dmrg iqdmrg dmrgj1j2

dmrg: dmrg.o $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) dmrg.o -o dmrg $(LIBFLAGS)

dmrg-g: mkdebugdir .debug_objs/dmrg.o $(ITENSOR_GLIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/dmrg.o -o dmrg-g $(LIBGFLAGS)

iqdmrg: iqdmrg.o $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) iqdmrg.o -o iqdmrg $(LIBFLAGS)

iqdmrg-g: mkdebugdir .debug_objs/iqdmrg.o $(ITENSOR_GLIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/iqdmrg.o -o iqdmrg-g $(LIBGFLAGS)

dmrgj1j2: dmrgj1j2.o $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) dmrgj1j2.o -o dmrgj1j2 $(LIBFLAGS)

dmrgj1j2-g: mkdebugdir .debug_objs/dmrgj1j2.o $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) .debug_objs/dmrgj1j2.o -o dmrgj1j2-g $(LIBFLAGS)

onesiteopt: onesiteopt.o $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) onesiteopt.o -o onesiteopt $(LIBFLAGS)

onesiteopt-g: mkdebugdir .debug_objs/onesiteopt.o $(ITENSOR_GLIBS) $(REL_TENSOR_HEADERS) 
	$(CCCOM) $(CCGFLAGS) .debug_objs/onesiteopt.o -o onesiteopt-g $(LIBGFLAGS)


mkdebugdir:
	mkdir -p .debug_objs

clean:
	rm -fr *.o .debug_objs dmrg dmrg-g iqdmrg iqdmrg-g onesiteopt onesiteopt-g \
	dmrgj1j2 dmrgj1j2-g
