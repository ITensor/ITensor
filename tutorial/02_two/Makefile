include ../../this_dir.mk
include ../../options.mk
################################################################

APP1=two
APP2=imag_tevol

LIBNAMES=itensor matrix utilities

#Define Dirs -----------
INCLUDEDIR=$(PREFIX)/include
INCLUDEFLAGS=. -I$(INCLUDEDIR) -I$(BOOST_DIR) $(BLAS_LAPACK_INCLUDEFLAGS)

#Mappings --------------
LOCAL_LIBFLAGS=$(patsubst %,-l%, $(LIBNAMES))
LOCAL_LIBGFLAGS=$(patsubst %,-l%-g, $(LIBNAMES))
LIBFILES=$(patsubst %,$(LIBDIR)/lib%.a, $(LIBNAMES))
LIBGFILES=$(patsubst %,$(LIBDIR)/lib%-g.a, $(LIBNAMES))
REL_TENSOR_HEADERS=$(patsubst %,$(INCLUDEDIR)/%, $(TENSOR_HEADERS))

#Define Flags ----------
CCFLAGS=$(CPPFLAGS) -I$(INCLUDEFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= -I$(INCLUDEFLAGS) -DDEBUG -DMATRIXBOUNDS -DITENSOR_USE_AT -DBOUNDS -g -Wall -ansi
LIBFLAGS=-L$(LIBDIR) $(LOCAL_LIBFLAGS) $(BLAS_LAPACK_LIBFLAGS)
LIBGFLAGS=-L$(LIBDIR) $(LOCAL_LIBGFLAGS) $(BLAS_LAPACK_LIBFLAGS)

#Rules ------------------

%.o: %.cc $(LIBFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(LIBFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: $(APP1)

debug: $(APP1)-g

$(APP1): $(APP1).o $(LIBFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) $(APP1).o -o $(APP1) $(LIBFLAGS)

$(APP1)-g: mkdebugdir .debug_objs/$(APP1).o $(LIBGFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/$(APP1).o -o $(APP1)-g $(LIBGFLAGS)

$(APP2): $(APP2).o $(LIBFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) $(APP2).o -o $(APP2) $(LIBFLAGS)

$(APP2)-g: mkdebugdir .debug_objs/$(APP2).o $(LIBGFILES) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/$(APP2).o -o $(APP2)-g $(LIBGFLAGS)

mkdebugdir:
	mkdir -p .debug_objs

clean:
	rm -fr *.o .debug_objs $(APP1) $(APP1)-g $(APP2) $(APP2)-g
