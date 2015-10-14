#
# Makefile for utilities library
#
#########################################

include ../this_dir.mk
include ../options.mk

#########################################

HEADERS=tinyformat.h print.h types.h indent.h cputime.h tarray1.h \
        flstring.h option.h input.h error.h minmax.h cplx_literal.h \
        autovector.h optional_ptr.h count.h stdx.h timers.h

OBJECTS= error.o ran1.o\
         cputime.o tarray1.o input.o option.o


CCFLAGS= -I. $(OPTIMIZATIONS) $(ITENSOR_INCLUDEFLAGS)
CCGFLAGS= -I. $(DEBUGFLAGS) $(ITENSOR_INCLUDEFLAGS)

GOBJECTS= $(patsubst %.o,.debug_objs/%.o, $(OBJECTS))

.SUFFIXES: .cc $(SUFFIXES)

.cc.o:	$*.cc
	@echo "Compiling utilities/$< with optimizations"
	$(eval COMMAND = $(CCCOM) -c $(CCFLAGS) -o $@ $<)
	@$(COMMAND) || (echo "Failure while executing command: $(COMMAND)" && exit 1)

.debug_objs/%.o: %.cc
	@echo "Compiling utilities/$< in debug mode"
	$(eval COMMAND = $(CCCOM) -c $(CCGFLAGS) -o $@ $<)
	@$(COMMAND) || (echo "Failure while executing command: $(COMMAND)" && exit 1)

install: headerscopy libutilities.a libutilities-g.a
	@cp libutilities.a libutilities-g.a $(ITENSOR_LIBDIR)
	@ranlib $(ITENSOR_LIBDIR)/libutilities.a
	@ranlib $(ITENSOR_LIBDIR)/libutilities-g.a

libutilities.a:	$(OBJECTS) $(HEADERS)
	@echo "Building static library libutilities.a"
	@ar r libutilities.a $(OBJECTS)
	@ranlib libutilities.a

libutilities-g.a: mkdebugdir $(GOBJECTS) $(HEADERS)
	@echo "Building static library libutilities-g.a"
	@ar r libutilities-g.a $(GOBJECTS)
	@ranlib libutilities-g.a

headerscopy:	
	@echo "Copying header files to $(ITENSOR_INCLUDEDIR)"
	@cp $(HEADERS) $(ITENSOR_INCLUDEDIR)

mkdebugdir:
	@mkdir -p .debug_objs

clean:	
	@echo "Cleaning build files in folder utilities/"
	@rm -f *.o libutilities.a libutilities-g.a
	@rm -fr .debug_objs
