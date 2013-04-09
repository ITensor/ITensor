#
# Makefile for Matrix library
#
#################################

include ../this_dir.mk
include ../options.mk

####################################################################

HEADERS=matrixref.h matrix.h precisio.h sparse.h bigmatrix.h davidson.h\
	storelink.h matrixref.ih matrix.ih conjugate_gradient.h sparseref.h\
    svd.h

OBJECTS=  matrix.o  utility.o  sparse.o  david.o sparseref.o\
	hpsortir.o  daxpy.o matrixref.o  storelink.o conjugate_gradient.o\
	 dgemm.o svd.o

SOURCES= matrix.cc utility.cc sparse.cc david.cc hpsortir.cc \
	matrixref.cc storelink.cc hpsortir.cc \
	conjugate_gradient.cc sparseref.cc\
	daxpy.cc svd.cc

GOBJECTS= $(patsubst %,g_objs/%, $(OBJECTS))

CCFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(OPTIMIZATIONS) -DPLATFORM_$(PLATFORM)
CCGFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS) -DPLATFORM_$(PLATFORM)

LIBFILE= libmatrix.a
LIBGFILE= libmatrix-g.a

.SUFFIXES: .cc $(SUFFIXES)

%.o: %.cc
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

g_objs/%.o: %.cc
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

install: libmatrix-g.a libmatrix.a
	cp libmatrix.a $(ITENSOR_LIBDIR)
	ranlib $(ITENSOR_LIBDIR)/libmatrix.a
	ranlib ./libmatrix.a
	cp libmatrix-g.a $(ITENSOR_LIBDIR)
	ranlib $(ITENSOR_LIBDIR)/libmatrix-g.a
	ranlib ./libmatrix-g.a
	cp *.h *.ih $(ITENSOR_INCLUDEDIR)

libmatrix.a: $(OBJECTS) $(HEADERS)
	ar rv libmatrix.a $(OBJECTS)
	touch libmatrix.a


libmatrix-g.a: mkgdir $(GOBJECTS) $(HEADERS)
	ar rv libmatrix-g.a $(GOBJECTS)
	touch libmatrix-g.a

mkgdir:
	mkdir -p g_objs

timeit:	timeit.o $(LIBFILE)
	$(CCCOM) $(CCFLAGS) -o timeit timeit.o $(ITENSOR_LIBFLAGS)
	touch timeit

bigtimeit:	bigtimeit.o $(LIBFILE)
	$(CCCOM) $(CCFLAGS) -o bigtimeit bigtimeit.o $(ITENSOR_LIBFLAGS)
	touch bigtimeit

testsref: testsref.o sparseref.o $(LIBFILE)
	$(CCCOM) $(CCGFLAGS) -o testsref testsref.o sparseref.o $(ITENSOR_LIBGFLAGS)
	touch testsref

testmem: g_objs/testmem.o $(LIBGFILE)
	$(CCCOM) $(CCGFLAGS) -o testmem g_objs/testmem.o $(ITENSOR_LIBGFLAGS)
	touch testmem

testmat:	g_objs/testmat.o $(LIBGFILE)
	$(CCCOM) $(CCGFLAGS) -o testmat g_objs/testmat.o $(ITENSOR_LIBGFLAGS)

clean:	
	rm -f *.o libmatrix.a libmatrix-g.a
	rm -fr g_objs
	rm -f test testsparse testsparse2 testsparse3 davidtest pbf testmat timeit
	rm -f subtest *.tar *.tar.gz

conjugate_gradient.o: matrix.h bigmatrix.h
sparseref.o: sparseref.h
storelink.o: storelink.h
matrixref.o: matrix.h matrixref.h storelink.h
matrix.o: matrix.h matrixref.h storelink.h
utility.o: matrix.h matrixref.h storelink.h
test.o: matrix.h matrixref.h storelink.h
testsparse.o: matrix.h sparse.h bigmatrix.h matrixref.h storelink.h
davidtest.o: matrix.h sparse.h bigmatrix.h davidson.h precisio.h\
	matrixref.h storelink.h
david.o: matrix.h sparse.h bigmatrix.h precisio.h matrixref.h storelink.h
sparse.o: matrix.h sparse.h bigmatrix.h matrixref.h storelink.h
svd.o: svd.h matrixref.h

g_objs/conjugate_gradient.o: matrix.h bigmatrix.h
g_objs/sparseref.o: sparseref.h
g_objs/storelink.o: storelink.h
g_objs/matrixref.o: matrix.h matrixref.h storelink.h
g_objs/matrix.o: matrix.h matrixref.h storelink.h
g_objs/utility.o: matrix.h matrixref.h storelink.h
g_objs/test.o: matrix.h matrixref.h storelink.h
g_objs/testsparse.o: matrix.h sparse.h bigmatrix.h matrixref.h storelink.h
g_objs/davidtest.o: matrix.h sparse.h bigmatrix.h davidson.h precisio.h\
	matrixref.h storelink.h
g_objs/david.o: matrix.h sparse.h bigmatrix.h precisio.h matrixref.h storelink.h
g_objs/sparse.o: matrix.h sparse.h bigmatrix.h matrixref.h storelink.h
g_objs/svd.o: svd.h matrixref.h
