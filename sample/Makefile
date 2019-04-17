include ../this_dir.mk
include ../options.mk

#Define Flags ----------

TENSOR_HEADERS=$(PREFIX)/itensor/all.h $(PREFIX)/itensor/mps/idmrg.h
CCFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(CPPFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS)
LIBFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBFLAGS)
LIBGFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBGFLAGS)

#Rules ------------------

%.o: %.cc $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: dmrg dmrg_table dmrgj1j2 exthubbard idmrg

debug: dmrg-g dmrg_table-g dmrgj1j2-g exthubbard-g idmrg-g

all: dmrg dmrg_table dmrgj1j2 exthubbard idmrg

dmrg: dmrg.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) dmrg.o -o dmrg $(LIBFLAGS)

dmrg-g: mkdebugdir .debug_objs/dmrg.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/dmrg.o -o dmrg-g $(LIBGFLAGS)

mixedspin: mixedspin.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) mixedspin.o -o mixedspin $(LIBFLAGS)

mixedspin-g: mkdebugdir .debug_objs/mixedspin.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/mixedspin.o -o mixedspin-g $(LIBGFLAGS)

dmrg_table: dmrg_table.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) dmrg_table.o -o dmrg_table $(LIBFLAGS)

dmrg_table-g: mkdebugdir .debug_objs/dmrg_table.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/dmrg_table.o -o dmrg_table-g $(LIBGFLAGS)

dmrgj1j2: dmrgj1j2.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) dmrgj1j2.o -o dmrgj1j2 $(LIBFLAGS)

dmrgj1j2-g: mkdebugdir .debug_objs/dmrgj1j2.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) .debug_objs/dmrgj1j2.o -o dmrgj1j2-g $(LIBGFLAGS)

exthubbard: exthubbard.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) exthubbard.o -o exthubbard $(LIBFLAGS)

exthubbard-g: mkdebugdir .debug_objs/exthubbard.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) .debug_objs/exthubbard.o -o exthubbard-g $(LIBGFLAGS)

idmrg: idmrg.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) idmrg.o -o idmrg $(LIBFLAGS)

idmrg-g: mkdebugdir .debug_objs/idmrg.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/idmrg.o -o idmrg-g $(LIBGFLAGS)

trg: trg.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) trg.o -o trg $(LIBFLAGS)

trg-g: mkdebugdir .debug_objs/trg.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/trg.o -o trg-g $(LIBGFLAGS)

mkdebugdir:
	mkdir -p .debug_objs

clean:
	@rm -fr *.o .debug_objs dmrg dmrg-g \
	dmrg_table dmrg_table-g dmrgj1j2 dmrgj1j2-g exthubbard exthubbard-g \
    idmrg idmrg-g mixedspin mixedspin-g trg trg-g
