LIBRARY_DIR=../..

APP1=two
APP2=imag_tevol

#################################################################
#################################################################
#################################################################
#################################################################

TENSOR_HEADERS=core.h

include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

#Mappings --------------
REL_TENSOR_HEADERS=$(patsubst %,$(ITENSOR_INCLUDEDIR)/%, $(TENSOR_HEADERS))

#Rules ------------------

%.o: %.cc $(REL_TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(REL_TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: $(APP1) $(APP2)
debug: $(APP1)-g $(APP2)-g

$(APP1): $(APP1).o $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $(APP1).o -o $(APP1) $(LIBFLAGS)

$(APP2): $(APP2).o $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $(APP2).o -o $(APP2) $(LIBFLAGS)

clean:
	rm -fr .debug_objs *.o $(APP1) $(APP1)-g $(APP2) $(APP2)-g

mkdebugdir:
	mkdir -p .debug_objs

