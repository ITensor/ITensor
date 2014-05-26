#
# Makefile for utilities library
#
#########################################

include ../this_dir.mk
include ../options.mk

#########################################

HEADERS=cppversion.h array2.h intarray2.h hash.h indent.h cputime.h tarray1.h \
        flstring.h option.h input.h error.h minmax.h safebool.h

OBJECTS= intarray2.o error.o ran1.o\
         hash.o cputime.o tarray1.o input.o option.o

SOURCES= array2.cc error.cc ran1.cc hash.cc\
         cputime.cc tarray1.cc testarray.cc intarray2.cc \
         input.cc option.cc

CCFLAGS= -I. $(OPTIMIZATIONS) $(ITENSOR_INCLUDEFLAGS)
CCGFLAGS= -I. $(DEBUGFLAGS) $(ITENSOR_INCLUDEFLAGS)

.SUFFIXES: .cc $(SUFFIXES)

.cc.o:	$*.cc
	$(CCCOM) -c $(CCFLAGS) $<

install: headerscopy libutilities.a libutilities-g.a
	cp libutilities.a libutilities-g.a $(ITENSOR_LIBDIR)
	ranlib $(ITENSOR_LIBDIR)/libutilities.a
	ranlib $(ITENSOR_LIBDIR)/libutilities-g.a

libutilities.a:	$(OBJECTS) $(HEADERS)
	ar rv libutilities.a $(OBJECTS)
	ranlib libutilities.a

libutilities-g.a: $(SOURCES) $(HEADERS)
	$(CCCOM) $(CCGFLAGS) -c intarray2.cc
	$(CCCOM) $(CCGFLAGS) -c tarray1.cc
	$(CCCOM) $(CCGFLAGS) -c error.cc
	$(CCCOM) $(CCGFLAGS) -c ran1.cc
	$(CCCOM) $(CCGFLAGS) -c hash.cc
	$(CCCOM) $(CCGFLAGS) -c cputime.cc
	$(CCCOM) $(CCGFLAGS) -c option.cc
	$(CCCOM) $(CCGFLAGS) -c input.cc
	ar rv libutilities-g.a $(OBJECTS)
	ranlib libutilities-g.a

headerscopy:	
	cp $(HEADERS) $(ITENSOR_INCLUDEDIR)

clean:	
	rm -f *.o libutilities.a libutilities-g.a
