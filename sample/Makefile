include ../options.mk
################################################################

TENSOR_HEADERS=tensor.h mps.h algorithm.h iq.h hams.h
LIBNAMES=matrix utilities itensor

#################################################################

#Define Dirs -----------
INCLUDEDIR=$(PREFIX)/include
INCLUDEFLAGS=. -I$(INCLUDEDIR) -I$(BOOST_DIR) $(BLAS_LAPACK_INCLUDEFLAGS)

#Mappings --------------
LOCAL_LIBFLAGS=$(patsubst %,-l%, $(LIBNAMES))
#LOCAL_LIBFLAGS=$(patsubst %,-l%-g, $(LIBNAMES))
LIBFILES=$(patsubst %,$(LIBDIR)/lib%.a, $(LIBNAMES))
#LIBFILES=$(patsubst %,$(LIBDIR)/lib%-g.a, $(LIBNAMES))
REL_TENSOR_HEADERS=$(patsubst %,$(INCLUDEDIR)/%, $(TENSOR_HEADERS))

#Define Flags ----------
CCFLAGS=$(CPPFLAGS) -I$(INCLUDEFLAGS) -O2  -DNDEBUG -DBOOST_DISABLE_ASSERTS
#CCFLAGS= -I$(INCLUDEFLAGS) -DMATRIXBOUNDS -DITENSOR_USE_AT -DBOUNDS -g -Wall -ansi
LIBFLAGS=-L$(LIBDIR) $(LOCAL_LIBFLAGS) $(BLAS_LAPACK_LIBFLAGS)

#Rules ------------------

%.o: %.cc $(REL_TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

#Targets -----------------

all: dmrg iqdmrg

dmrg: dmrg.o $(LIBFILES)
	$(CCCOM) $(CCFLAGS) dmrg.o -o dmrg $(LIBFLAGS)

iqdmrg: iqdmrg.o $(LIBFILES)
	$(CCCOM) $(CCFLAGS) iqdmrg.o -o iqdmrg $(LIBFLAGS)

clean:
	rm -f *.o dmrg iqdmrg
