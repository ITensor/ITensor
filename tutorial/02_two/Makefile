include ../../this_dir.mk
include ../../options.mk
################################################################

APP1=two
APP2=imag_tevol

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

build: $(APP1) $(APP2)

debug: $(APP1)-g $(APP2)-g

$(APP1): $(APP1).o $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) $(APP1).o -o $(APP1) $(LIBFLAGS)

$(APP1)-g: mkdebugdir .debug_objs/$(APP1).o $(ITENSOR_GLIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/$(APP1).o -o $(APP1)-g $(LIBGFLAGS)

$(APP2): $(APP2).o $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) $(APP2).o -o $(APP2) $(LIBFLAGS)

$(APP2)-g: mkdebugdir .debug_objs/$(APP2).o $(ITENSOR_GLIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/$(APP2).o -o $(APP2)-g $(LIBGFLAGS)

mkdebugdir:
	mkdir -p .debug_objs

clean:
	rm -fr *.o .debug_objs $(APP1) $(APP1)-g $(APP2) $(APP2)-g
