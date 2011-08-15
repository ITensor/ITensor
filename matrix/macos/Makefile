#
# Makefile for Matrix library
#
#################################

include ../../this_dir.mk
include ../../options.mk

####################################################################

HEADERS=matrixref.h matrix.h precisio.h sparse.h bigmatrix.h davidson.h\
	storelink.h matrixref.ih matrix.ih conjugate_gradient.h sparseref.h

OBJECTS=  matrix.o  utility.o  sparse.o  david.o sparseref.o\
	hpsortir.o  daxpy.o matrixref.o  storelink.o conjugate_gradient.o\
	 dgemm.o 

OOBJECTS=  matrix.o-o  utility.o-o  sparse.o-o  david.o-o sparseref.o-o\
	daxpy.o-o hpsortir.o-o conjugate_gradient.o-o  matrixref.o-o  storelink.o-o \
	dgemm.o-o

SOURCES= matrix.cc utility.cc sparse.cc david.cc hpsortir.cc \
	matrixref.cc storelink.cc hpsortir.cc \
	conjugate_gradient.cc sparseref.cc\
	daxpy.cc 

GOBJECTS= $(patsubst %,.g_objs/%, $(OBJECTS))
PGOBJECTS= $(patsubst %,.pg_objs/%, $(OBJECTS))

INCLUDEFLAGS=. -I$(INCLUDEDIR) $(BLAS_LAPACK_INCLUDEFLAGS)

CCFLAGS= -I$(INCLUDEFLAGS) -O2 -DNDEBUG
CCGFLAGS= -I$(INCLUDEFLAGS) -DMATRIXBOUNDS -DBOUNDS -g -O0
CCPGFLAGS= -I$(INCLUDEFLAGS) -pg -O

LIBFLAGS= -L$(LIBDIR) -lmatrix -lutilities $(BLAS_LAPACK_LIBFLAGS)
LIBGFLAGS= -L$(LIBDIR) -lmatrix-g -lutilities-g $(BLAS_LAPACK_LIBFLAGS)

LIBFILE= libmatrix.a
LIBGFILE= libmatrix-g.a

.SUFFIXES: .cc $(SUFFIXES) .o-g .o-pg .o-o

%.o: %.cc
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.g_objs/%.o: %.cc
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

.pg_objs/%.o: %.cc
	$(CCCOM) -c $(CCPGFLAGS) -o $@ $<

install: libmatrix-g.a libmatrix.a
	cp libmatrix.a $(LIBDIR)
	ranlib $(LIBDIR)/libmatrix.a
	ranlib ./libmatrix.a
	cp libmatrix-g.a $(LIBDIR)
	ranlib $(LIBDIR)/libmatrix-g.a
	ranlib ./libmatrix-g.a
	cp *.h *.ih $(INCLUDEDIR)

libmatrix.a: $(OBJECTS) $(HEADERS)
	ar rv libmatrix.a $(OBJECTS)
	touch libmatrix.a


libmatrix-g.a: mkgdir $(GOBJECTS) $(HEADERS)
	ar rv libmatrix-g.a $(GOBJECTS)
	touch libmatrix-g.a

mkgdir:
	mkdir -p .g_objs

libmatrix-pg.a: mkpgdir $(PGOBJECTS) $(HEADERS)
	ar rv libmatrix-pg.a $(PGOBJECTS)
	touch libmatrix-pg.a

mkpgdir:
	mkdir -p .pg_objs

timeit:	timeit.o-o $(LIBFILE)
	mv timeit.o-o timeit.o
	$(CCCOM) $(CCFLAGS) -o timeit timeit.o $(LIBFLAGS)
	mv timeit.o timeit.o-o
	touch timeit

bigtimeit:	bigtimeit.o-o $(LIBFILE)
	mv bigtimeit.o-o bigtimeit.o
	$(CCCOM) $(CCFLAGS) -o bigtimeit bigtimeit.o $(LIBFLAGS)
	mv bigtimeit.o bigtimeit.o-o
	touch bigtimeit

testsref:	testsref.o-o sparseref.o-o $(LIBGFILE)
	mv testsref.o-o testsref.o
	$(CCCOM) $(CCGFLAGS) -o testsref testsref.o sparseref.o-o $(LIBGFLAGS)
	mv testsref.o testsref.o-o
	touch testsref

testmem:	testmem.o-g $(LIBGFILE)
	mv testmem.o-g testmem.o
	$(CCCOM) $(CCGFLAGS) -o testmem testmem.o $(LIBGFLAGS)
	mv testmem.o testmem.o-g
	touch testmem

testmat:	testmat.o-g $(LIBGFILE)
	mv testmat.o-g testmat.o
#	$(CCCOM) $(CCGFLAGS) -o testmat testmat.o cblas_dgemm.o $(LIBGFLAGS)
	$(CCCOM) $(CCGFLAGS) -o testmat testmat.o $(LIBGFLAGS)
#-L. -lmatrix-g -lutilities-g -framework Accelerate
#-L. -lmatrix-g -L./lib -lutilities -lcblas -latlas -lm

clean:	
	rm -f *.o *.o-g *.o-o *.o-pg libmatrix.a libmatrix-g.a libmatrix-pg.a 
	rm -fr .g_objs .pg_objs
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

.g_objs/conjugate_gradient.o: matrix.h bigmatrix.h
.g_objs/sparseref.o: sparseref.h
.g_objs/storelink.o: storelink.h
.g_objs/matrixref.o: matrix.h matrixref.h storelink.h
.g_objs/matrix.o: matrix.h matrixref.h storelink.h
.g_objs/utility.o: matrix.h matrixref.h storelink.h
.g_objs/test.o: matrix.h matrixref.h storelink.h
.g_objs/testsparse.o: matrix.h sparse.h bigmatrix.h matrixref.h storelink.h
.g_objs/davidtest.o: matrix.h sparse.h bigmatrix.h davidson.h precisio.h\
	matrixref.h storelink.h
.g_objs/david.o: matrix.h sparse.h bigmatrix.h precisio.h matrixref.h storelink.h
.g_objs/sparse.o: matrix.h sparse.h bigmatrix.h matrixref.h storelink.h

.pg_objs/conjugate_gradient.o: matrix.h bigmatrix.h
.pg_objs/sparseref.o: sparseref.h
.pg_objs/storelink.o: storelink.h
.pg_objs/matrixref.o: matrix.h matrixref.h storelink.h
.pg_objs/matrix.o: matrix.h matrixref.h storelink.h
.pg_objs/utility.o: matrix.h matrixref.h storelink.h
.pg_objs/test.o: matrix.h matrixref.h storelink.h
.pg_objs/testsparse.o: matrix.h sparse.h bigmatrix.h matrixref.h storelink.h
.pg_objs/davidtest.o: matrix.h sparse.h bigmatrix.h davidson.h precisio.h\
	matrixref.h storelink.h
.pg_objs/david.o: matrix.h sparse.h bigmatrix.h precisio.h matrixref.h storelink.h
.pg_objs/sparse.o: matrix.h sparse.h bigmatrix.h matrixref.h storelink.h
	
#dependencies

conjugate_gradient.o-o: matrix.h bigmatrix.h

conjugate_gradient.o-g: matrix.h bigmatrix.h

conjugate_gradient.o-pg: matrix.h bigmatrix.h

sparseref.o-o: sparseref.h

sparseref.o-g: sparseref.h

sparseref.o-pg: sparseref.h

testsref.o-g: sparseref.h

storelink.o-o: storelink.h

storelink.o-g: storelink.h

storelink.o-pg: storelink.h

matrixref.h: matrixref.ih
	touch matrixref.h

matrix.h: matrix.ih
	touch matrix.h

matrixref.o-o: matrix.h matrixref.h storelink.h

matrixref.o-g: matrix.h matrixref.h storelink.h

matrixref.o-pg: matrix.h matrixref.h storelink.h

matrix.o-o: matrix.h matrixref.h storelink.h

matrix.o-g: matrix.h matrixref.h storelink.h

matrix.o-pg: matrix.h matrixref.h storelink.h

utility.o-o: matrix.h matrixref.h storelink.h

utility.o-g: matrix.h matrixref.h storelink.h

utility.o-pg: matrix.h matrixref.h storelink.h

test.o-o: matrix.h matrixref.h storelink.h

test.o-g: matrix.h matrixref.h storelink.h

test.o-pg: matrix.h matrixref.h storelink.h

testsparse.o-o: matrix.h sparse.h bigmatrix.h matrixref.h storelink.h

testsparse.o-g: matrix.h sparse.h bigmatrix.h matrixref.h storelink.h

testsparse.o-pg: matrix.h sparse.h bigmatrix.h matrixref.h storelink.h

testmat.o-g: matrix.h sparse.h bigmatrix.h matrixref.h storelink.h

davidtest.o-o: matrix.h sparse.h bigmatrix.h davidson.h precisio.h\
	matrixref.h storelink.h

davidtest.o-g: matrix.h sparse.h bigmatrix.h davidson.h precisio.h\
	matrixref.h storelink.h

davidtest.o-pg: matrix.h sparse.h bigmatrix.h davidson.h precisio.h\
	matrixref.h storelink.h

david.o-o: matrix.h sparse.h bigmatrix.h precisio.h matrixref.h storelink.h

david.o-g: matrix.h sparse.h bigmatrix.h precisio.h matrixref.h storelink.h

david.o-pg: matrix.h sparse.h bigmatrix.h precisio.h matrixref.h\
	storelink.h

sparse.o-o: matrix.h sparse.h bigmatrix.h matrixref.h storelink.h

sparse.o-g: matrix.h sparse.h bigmatrix.h matrixref.h storelink.h

sparse.o-pg: matrix.h sparse.h bigmatrix.h matrixref.h storelink.h
