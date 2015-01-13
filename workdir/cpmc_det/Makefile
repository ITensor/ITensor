# 1. Rename this file to just 'Makefile' if it's not called that already. 
#    Put it in the same folder as your 'driver' code (the code containing
#    the 'main' function).

# 2. Edit LIBRARY_DIR to point at the location of your ITensor Library
#    source folder (this is the folder that has options.mk in it)
LIBRARY_DIR=$(HOME)/ITensor

# 3. If your 'main' function is in a file called 'myappname.cc', then
#    set APP to 'myappname'. Running 'make' will compile the app.
#    Running 'make debug' will make a program called 'myappname-g'
#    which includes debugging symbols and can be used in gdb (Gnu debugger);
APP=run_cpmc

# 4. Add any headers your program depends on here. The make program
#    will auto-detect if these headers have changed and recompile your app.
HEADERS=cpmc.h

# 5. For any additional .cc files making up your project,
#    add their full filenames here.
CCFILES=$(APP).cc cpmc.cc

#################################################################
#################################################################
#################################################################
#################################################################

TENSOR_HEADERS=core.h

include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

#Mappings --------------
REL_TENSOR_HEADERS=$(patsubst %,$(ITENSOR_INCLUDEDIR)/%, $(TENSOR_HEADERS))
OBJECTS=$(patsubst %.cc,%.o, $(CCFILES))
GOBJECTS=$(patsubst %,.debug_objs/%, $(OBJECTS))

#Rules ------------------

%.o: %.cc $(HEADERS) $(REL_TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(HEADERS) $(REL_TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: $(APP)
debug: $(APP)-g

$(APP): $(OBJECTS) $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $(OBJECTS) -o $(APP) $(LIBFLAGS)

$(APP)-g: mkdebugdir $(GOBJECTS) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(APP)-g $(LIBGFLAGS)

clean:
	rm -fr .debug_objs *.o $(APP) $(APP)-g

mkdebugdir:
	mkdir -p .debug_objs

