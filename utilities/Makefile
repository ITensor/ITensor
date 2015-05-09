#
# Makefile for utilities library
#
#########################################

include ../this_dir.mk
include ../options.mk

#########################################

HEADERS=tinyformat.h print.h types.h indent.h cputime.h \
        flstring.h args.h input.h error.h  cplx_literal.h \
        autovector.h optional_ptr.h count.h

OBJECTS = error.o 
OBJECTS += cputime.o 
OBJECTS += input.o 
OBJECTS += args.o

CCFLAGS= -I. $(OPTIMIZATIONS) $(ITENSOR_INCLUDEFLAGS) -DPLATFORM_$(PLATFORM)
CCGFLAGS= -I. $(DEBUGFLAGS) $(ITENSOR_INCLUDEFLAGS) -DPLATFORM_$(PLATFORM)

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
