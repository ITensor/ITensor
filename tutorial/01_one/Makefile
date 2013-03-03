include ../../this_dir.mk
include ../../options.mk
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

build: one

debug: one-g

one: one.o $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) one.o -o one $(LIBFLAGS)

one-g: mkdebugdir .debug_objs/one.o $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/one.o -o one-g $(LIBGFLAGS)

mkdebugdir:
	mkdir -p .debug_objs

clean:
	rm -fr *.o .debug_objs one one-g
