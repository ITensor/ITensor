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

SOURCES= matrix.cc utility.cc sparse.cc david.cc hpsortir.cc \
	matrixref.cc storelink.cc hpsortir.cc \
	conjugate_gradient.cc sparseref.cc\
	daxpy.cc 

GOBJECTS= $(patsubst %,.g_objs/%, $(OBJECTS))
PGOBJECTS= $(patsubst %,.pg_objs/%, $(OBJECTS))

INCLUDEFLAGS=. -I$(INCLUDEDIR) $(BLAS_LAPACK_INCLUDEFLAGS)
CCFLAGS= -I$(INCLUDEFLAGS) $(OPTIMIZATION) 
CCGFLAGS= -I$(INCLUDEFLAGS) -DMATRIXBOUNDS -DBOUNDS -g -O0
CCPGFLAGS= -I$(INCLUDEFLAGS) -pg -O

LIBFLAGS= -L$(LIBDIR) -L$(BLAS_LAPACK_LIBDIR) -lmatrix -lutilities $(BLAS_LAPACK_FLAGS)
LIBGFLAGS= -L$(LIBDIR) -L$(BLAS_LAPACK_LIBDIR) -lmatrix-g -lutilities-g $(BLAS_LAPACK_FLAGS)

LIBFILE=libmatrix.a
LIBGFILE=libmatrix-g.a

.SUFFIXES: .cc $(SUFFIXES) .o-g .o-pg .o-o

.cc.o:	$*.cc	Makefile
	$(CCCOM) -c $(CCFLAGS) $<
	touch $*.o

.cc.s:	$*.cc	Makefile
	$(CCCOM) -c -S -O1 -I$(INCLUDEDIR) $<
#	unop.p $*.o

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

daxpy.o-o: daxpy.cc
	$(CCCOM) $(OPTIMIZATION) -c daxpy.cc
	mv daxpy.o daxpy.o-o

daxpy.o-g: daxpy.cc
	$(CCCOM) $(OPTIMIZATION) -c daxpy.cc
	mv daxpy.o daxpy.o-g

#dgemm.o-o: dgemm.cc
#	$(CCCOM) -O1 -I$(INCLUDEDIR) -S -c dgemm.cc
#	unop.p dgemm.s
#	$(CCCOM) -O1 -I$(INCLUDEDIR) -c dgemm.s
#	mv dgemm.o dgemm.o-o

#utility.o-o: utility.cc
#	$(CCCOM) -O1 -I$(INCLUDEDIR) -S -c utility.cc
#	unop.p utility.s
#	$(CCCOM) -O1 -I$(INCLUDEDIR) -c utility.s
#	mv utility.o utility.o-o

#dgemm.o-g: dgemm.cc
#	$(CCCOM) -g -O1 -I$(INCLUDEDIR) -S -c dgemm.cc
#	unop.p dgemm.s
#	$(CCCOM) -g -O1 -I$(INCLUDEDIR) -c dgemm.s
#	mv dgemm.o dgemm.o-g

donum:	donum.o-g $(LIBGFILE)
	mv donum.o-g donum.o
	$(CCCOM) $(CCGFLAGS) -o donum donum.o $(LIBGFLAGS)
	mv donum.o donum.o-g
	touch donum

shank:	shank.o-g $(LIBGFILE)
	mv shank.o-g shank.o
	$(CCCOM) $(CCGFLAGS) -o shank shank.o $(LIBGFLAGS)
	mv shank.o shank.o-g
	touch shank

lorentzfit:	lorentzfit.o-g $(LIBGFILE)
	mv lorentzfit.o-g lorentzfit.o
	$(CCCOM) $(CCGFLAGS) -o lorentzfit lorentzfit.o $(LIBGFLAGS)
	mv lorentzfit.o lorentzfit.o-g
	touch lorentzfit

linfit2:	linfit2.o-g $(LIBGFILE)
	mv linfit2.o-g linfit2.o
	$(CCCOM) $(CCGFLAGS) -o linfit2 linfit2.o $(LIBGFLAGS)
	mv linfit2.o linfit2.o-g
	touch linfit2

dendiag:	dendiag.o-g $(LIBGFILE)
	mv dendiag.o-g dendiag.o
	$(CCCOM) $(CCGFLAGS) -o dendiag dendiag.o $(LIBGFLAGS)
	mv dendiag.o dendiag.o-g
	touch dendiag

doft:	doft.o-g $(LIBGFILE)
	mv doft.o-g doft.o
	$(CCCOM) $(CCGFLAGS) -o doft doft.o $(LIBGFLAGS)
	mv doft.o doft.o-g
	touch doft

fun:	fun.o-g $(LIBGFILE)
	mv fun.o-g fun.o
	$(CCCOM) $(CCGFLAGS) -o fun fun.o $(LIBGFLAGS)
	mv fun.o fun.o-g
	touch fun

emat:	emat.o-g $(LIBGFILE)
	mv emat.o-g emat.o
	$(CCCOM) $(CCGFLAGS) -o emat emat.o $(LIBGFLAGS)
	mv emat.o emat.o-g
	touch emat

imat:	imat.o-g $(LIBGFILE)
	mv imat.o-g imat.o
	$(CCCOM) $(CCGFLAGS) -o imat imat.o $(LIBGFLAGS)
	mv imat.o imat.o-g
	touch imat

timeit:	timeit.o-o $(LIBFILE)
	mv timeit.o-o timeit.o
	$(CCCOM) $(CCFLAGS) -o timeit timeit.o $(LIBFLAGS)
	mv timeit.o timeit.o-o
	touch timeit

timeit2:	timeit2.o-o $(LIBFILE)
	mv timeit2.o-o timeit2.o
	$(CCCOM) $(CCFLAGS) -o timeit2 timeit2.o $(LIBFLAGS)
	mv timeit2.o timeit2.o-o
	touch timeit2

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

#	$(CCCOM) $(CCGFLAGS) -o testmat testmat.o lib/tmp/cblas_dgemm.o $(LIBGFLAGS)
	mv testmat.o testmat.o-g
	touch testmat

testeig:	testeig.o-g $(LIBFILE)
	mv testeig.o-g testeig.o
	$(CCCOM) $(CCFLAGS) -o testeig testeig.o $(LIBFLAGS)
	mv testeig.o testeig.o-g
	touch testeig

testeign:	testeig.o-g $(LIBFILE) ut2.o-o
	ln -f testeig.o-g testeig.o
	ln -f ut2.o-o ut2.o
	$(CCCOM) $(CCFLAGS) -g -o testeign testeig.o ut2.o $(LIBFLAGS)
	touch testeign

convolve:	convolve.o-g $(LIBGFILE)
	mv convolve.o-g convolve.o
	$(CCCOM) $(CCGFLAGS) -o convolve convolve.o $(LIBGFLAGS)
	mv convolve.o convolve.o-g
	touch convolve

testquick:	testquick.o $(LIBFILE)
	$(CCCOM) $(CCFLAGS) -o testquick testquick.o $(LIBFLAGS)

matrixref.trans:	matrixref.cc matrixref.h matrixref.ih
	$(CCCOM) -F -I$(INCLUDEDIR) matrixref.cc | c++filt > matrixref.trans

storelink.trans:	storelink.cc storelink.h
	CC -F -I$(INCLUDEDIR) storelink.cc | c++filt > storelink.trans

subtest: subtest.o $(LIBFILE)
	$(CCCOM) $(CCFLAGS) -o subtest subtest.o $(LIBFLAGS)

testsparse:	testsparse.o $(LIBFILE)
	$(CCCOM) $(CCFLAGS) -o testsparse testsparse.o $(LIBFLAGS) 

testsparse2:	testsparse2.o $(LIBFILE)
	$(CCCOM) $(CCFLAGS) -o testsparse2 testsparse2.o $(LIBFLAGS) 

testsparse3:	testsparse3.o $(LIBFILE)
	$(CCCOM) $(CCFLAGS) -o testsparse3 testsparse3.o $(LIBFLAGS) 

sparsemult:	sparsemult.o $(LIBFILE)
	$(CCCOM) $(CCFLAGS) -o sparsemult sparsemult.o $(LIBFLAGS) 

davidtest: davidtest.o-g $(LIBGFILE)
	mv davidtest.o-g davidtest.o
	$(CCCOM) $(CCGFLAGS) -o davidtest davidtest.o $(LIBGFLAGS)
	mv davidtest.o davidtest.o-g
	touch davidtest

bench:	bench.o-o $(LIBFILE)
	mv bench.o-o bench.o
	$(CCCOM) $(CCFLAGS) -o bench bench.o $(LIBFLAGS)
	mv bench.o bench.o-o
	touch bench

pbf:	pbf.o
	$(CCCOM) $(CCFLAGS) -o pbf pbf.o $(LIBFLAGS)

clean:	
	rm -f *.o *.o-g *.o-o *.o-pg libmatrix.a libmatrix-g.a libmatrix-pg.a 
	rm -f test testsparse testsparse2 testsparse3 davidtest pbf timeit testmat
	rm -f subtest *.tar *.tar.gz


c90:	tar
	put_c90 matrix.tar
	put_c90 include.tar

targz:	matrix.tar.gz include.tar.gz

tar:	matrix.tar include.tar

matrix.tar: $(HEADERS) $(SOURCES) 
	cp Makefile Makefile.new
	tar cvf ../matrix.tar $(HEADERS) $(SOURCES) Makefile.new
	mv ../matrix.tar .

matrix.tar.gz: matrix.tar
	gzip -f matrix.tar

include.tar: *.h Makefile
	tar cvf ../include.tar -C $(INCLUDEDIR) .
	mv ../include.tar .

include.tar.gz: include.tar
	gzip -f include.tar

sdsc: matrix.tar #include.tar
	putc90 matrix.tar matrix/matrix.tar
#	putc90 include.tar matrix/include.tar

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
