include ../this_dir.mk
include ../options.mk

#Define Flags ----------

TENSOR_HEADERS=$(PREFIX)/itensor/all.h
CCFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(CPPFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS)
LIBFLAGS=-L'$(ITENSOR_LIBDIR)' $(ITENSOR_LIBFLAGS)
LIBGFLAGS=-L'$(ITENSOR_LIBDIR)' $(ITENSOR_LIBGFLAGS)

#Rules ------------------

%.o: %.cc $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: dmrg dmrg_table dmrgj1j2 exthubbard trg ctmrg mixedspin hubbard_2d hubbard_2d_conserve_momentum

debug: dmrg-g dmrg_table-g dmrgj1j2-g exthubbard-g trg-g ctmrg-g mixedspin-g hubbard_2d-g hubbard_2d_conserve_momentum-g

all: dmrg dmrg_table dmrgj1j2 exthubbard trg ctmrg mixedspin hubbard_2d hubbard_2d_conserve_momentum

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

trg: trg.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) trg.o -o trg $(LIBFLAGS)

trg-g: mkdebugdir .debug_objs/trg.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/trg.o -o trg-g $(LIBGFLAGS)

ctmrg: ctmrg.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) ctmrg.o -o ctmrg $(LIBFLAGS)

ctmrg-g: mkdebugdir .debug_objs/trg.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/ctmrg.o -o ctmrg-g $(LIBGFLAGS)

hubbard_2d: hubbard_2d.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) hubbard_2d.o -o hubbard_2d $(LIBFLAGS)

hubbard_2d-g: mkdebugdir .debug_objs/hubbard_2d.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/hubbard_2d.o -o hubbard_2d-g $(LIBGFLAGS)

hubbard_2d_conserve_momentum: hubbard_2d_conserve_momentum.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) hubbard_2d_conserve_momentum.o -o hubbard_2d_conserve_momentum $(LIBFLAGS)

hubbard_2d_conserve_momentum-g: mkdebugdir .debug_objs/hubbard_2d_conserve_momentum.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/hubbard_2d_conserve_momentum.o -o hubbard_2d_conserve_momentum-g $(LIBGFLAGS)

mkdebugdir:
	mkdir -p .debug_objs

clean:
	@rm -fr *.o .debug_objs dmrg dmrg-g \
	dmrg_table dmrg_table-g dmrgj1j2 dmrgj1j2-g exthubbard exthubbard-g \
    mixedspin mixedspin-g trg trg-g ctmrg ctmrg-g hubbard_2d hubbard_2d-g hubbard_2d_conserve_momentum hubbard_2d_conserve_momentum-g
