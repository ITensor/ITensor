include ../../this_dir.mk
include ../../options.mk
################################################################

APP=dmrg

################################################################

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

build: $(APP)

debug: $(APP)-g

$(APP): $(APP).o $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) $(APP).o -o $(APP) $(LIBFLAGS)

$(APP)-g: mkdebugdir .debug_objs/$(APP).o $(ITENSOR_GLIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/$(APP).o -o $(APP)-g $(LIBGFLAGS)

mkdebugdir:
	mkdir -p .debug_objs

clean:
	rm -fr *.o .debug_objs $(APP) $(APP)-g
