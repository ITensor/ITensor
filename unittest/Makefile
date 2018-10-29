include ../this_dir.mk
include ../options.mk
################################################################

HEADR=$(ITENSOR_INCLUDEDIR)/itensor

HEADERS=test.h
SOURCES=test.cc
ifdef test
SOURCES+= $(test)_test.cc
else
SOURCES+= util_test.cc
SOURCES+= algorithm_test.cc
SOURCES+= real_test.cc
SOURCES+= args_test.cc
SOURCES+= matrix_test.cc
SOURCES+= tensor_test.cc
#SOURCES+= contract_test.cc
SOURCES+= sparse_contract_test.cc
SOURCES+= index_test.cc
SOURCES+= indexset_test.cc
SOURCES+= itensor_test.cc
SOURCES+= qn_test.cc
SOURCES+= iqindex_test.cc
SOURCES+= iqtensor_test.cc
SOURCES+= decomp_test.cc
SOURCES+= mps_test.cc
SOURCES+= mpo_test.cc
SOURCES+= autompo_test.cc
SOURCES+= iterativesolvers_test.cc
SOURCES+= regression_test.cc
#SOURCES+= spectrum_test.cc
#SOURCES+= webpage_test.cc
SOURCES+= localop_test.cc
SOURCES+= siteset_test.cc
#SOURCES+= bondgate_test.cc
endif

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
	@echo "Compiling $< with optimizations"
	$(eval COMMAND = $(CCCOM) -c $(CCFLAGS) -o $@ $<)
	@$(COMMAND) || (echo "Failure while executing command: $(COMMAND)" && exit 1)

.debug_objs/%.o: %.cc $(ITENSOR_GLIBS)
	@echo "Compiling $< in debug mode"
	$(eval COMMAND = $(CCCOM) -c $(CCGFLAGS) -o $@ $<)
	@$(COMMAND) || (echo "Failure while executing command: $(COMMAND)" && exit 1)

#Targets -----------------
#@$(CCCOM) $(CCGFLAGS) .debug_objs/test.o .debug_objs/%_test.o -o test-g $(LIBGFLAGS)

test_%: mkdebugdir .debug_objs/test.o $(ITENSOR_GLIBS) .debug_objs/%_test.o
	@echo 
	@echo Running test index...
	@echo 
	@./test-g

run: test-g
	@echo 
	@echo Running all tests...
	@echo 
	@./test-g

debug: test-g

test-g: mkdebugdir .debug_objs/test.o $(ITENSOR_GLIBS) $(GOBJECTS)
	@$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o test-g $(LIBGFLAGS)


mkdebugdir:
	@mkdir -p .debug_objs

clean:
	@rm -fr *.o .debug_objs test test-g


LIBHEADERS=$(HEADR)/util/infarray.h
util_test.cc: $(LIBHEADERS)
.debug/util_test.cc: $(LIBHEADERS)

LIBHEADERS=$(HEADR)/matrix/vec.h $(HEADR)/matrix/mat.h
matrix_test.o: $(LIBHEADERS)
.debug/matrix_test.o: $(LIBHEADERS)

LIBHEADERS=$(ITENSOR_INCLUDEDIR)/itensor/real.h
real_test.o: $(LIBHEADERS)
.debug/real_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/itensor/index.h
index_test.o: $(LIBHEADERS)
.debug_objs/index_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/itensor/itensor.h
itensor_test.o: $(LIBHEADERS)
.debug_objs/itensor_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/itensor/mps/mps.h
mps_test.o: $(LIBHEADERS)
.debug_objs/mps_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/itensor/mps/mpo.h
mpo_test.o: $(LIBHEADERS)
.debug_objs/mpo_test.o: $(LIBHEADERS)

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/itensor/decomp.h
decomp_test.o: $(LIBHEADERS)
.debug_objs/decomp_test.o: $(LIBHEADERS)

regression_test.o: $(LIBHEADERS)
.debug_objs/regression_test.o: $(LIBHEADERS)

iterativesolvers_test.o: $(ITENSOR_INCLUDEDIR)/itensor/iterativesolvers.h
.debug_objs/iterativesolvers_test.o: $(ITENSOR_INCLUDEDIR)/itensor/iterativesolvers.h

webpage_test.o: $(ITENSOR_INCLUDEDIR)/itensor/iqtensor.h
.debug_objs/webpage_test.o: $(ITENSOR_INCLUDEDIR)/itensor/iqtensor.h

LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/itensor/mps/localop.h
LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/itensor/mps/localmpo.h
LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/itensor/mps/localmpo_mps.h
LIBHEADERS+= $(ITENSOR_INCLUDEDIR)/itensor/mps/localmposet.h
localop_test.o: $(LIBHEADERS)
.debug_objs/localop_test.o: $(LIBHEADERS)

