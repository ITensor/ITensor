#
# Makefile for utilities library
#
#########################################

include ../this_dir.mk
include ../options.mk

#########################################

HEADERS=tinyformat.h print.h cppversion.h indent.h cputime.h tarray1.h \
        flstring.h option.h input.h error.h minmax.h safebool.h

OBJECTS= error.o ran1.o\
         cputime.o tarray1.o input.o option.o


CCFLAGS= -I. $(OPTIMIZATIONS) $(ITENSOR_INCLUDEFLAGS)
CCGFLAGS= -I. $(DEBUGFLAGS) $(ITENSOR_INCLUDEFLAGS)

GOBJECTS= $(patsubst %.o,.debug_objs/%.o, $(OBJECTS))

.SUFFIXES: .cc $(SUFFIXES)

.cc.o:	$*.cc
	$(CCCOM) -c $(CCFLAGS) $<

.debug_objs/%.o: %.cc
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

install: headerscopy libutilities.a libutilities-g.a
	cp libutilities.a libutilities-g.a $(ITENSOR_LIBDIR)
	ranlib $(ITENSOR_LIBDIR)/libutilities.a
	ranlib $(ITENSOR_LIBDIR)/libutilities-g.a

libutilities.a:	$(OBJECTS) $(HEADERS)
	ar rv libutilities.a $(OBJECTS)
	ranlib libutilities.a

libutilities-g.a: mkdebugdir $(GOBJECTS) $(HEADERS)
	ar rv libutilities-g.a $(GOBJECTS)
	ranlib libutilities-g.a

headerscopy:	
	cp $(HEADERS) $(ITENSOR_INCLUDEDIR)

mkdebugdir:
	mkdir -p .debug_objs

clean:	
	rm -f *.o libutilities.a libutilities-g.a
	rm -fr .debug_objs
