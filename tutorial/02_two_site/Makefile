LIBRARY_DIR=../..

APP1=two
APP2=imag_tevol

#################################################################
#################################################################
#################################################################
#################################################################

include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk


#Rules ------------------

%.o: %.cc
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc
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

