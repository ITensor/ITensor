// matrixref.cc -- Code for MatrixRef class

#include "matrix.h"
#include <math.h>
//#include <stdlib.h>
#include <iomanip>
#include <memory>
#include "indent.h"
#if defined(i386) || defined(__x86_64)
extern "C" void dcopy_(int*,Real*,int*,Real*,int*);
extern "C" void dscal_(int*,Real*,Real*,int*);
extern "C" void	daxpy_(int*,Real*,Real*,int*,Real*,int*);
extern "C" void dset_(int*,Real*,Real*,int*);
extern "C" void dvcal_(int*,Real*,Real*,int*,Real*,int*);
extern "C" Real ddot_(int*,Real*,int*,Real*,int*);
extern "C" Real dnrm2_(int*,Real*,int*);
extern "C" Real dsum_(int*,Real*,int*);
extern "C" void dgemt_(char*,int*,int*,Real*,Real*,int*,Real*,int*);
extern "C" void dgema_(char*,char*,int*,int*,Real*,Real*,int*,
			Real*,Real*,int*,Real*,int*);
extern "C" void dgemv_(char*,int*,int*,Real*,Real*,int*,Real*,int*,
				    Real*,Real*,int*);
extern "C" void dgemm_(char*,char*,int*,int*,int*,Real*,Real*,int*,
				Real*,int*,Real*,Real*,int*);
#else
void daxpy(int n, Real alpha, Real* x, int incx, Real* y, int incy);

#endif


// Print an error message and abort for debugging
void 
_merror(const char *s)
    {
    cerr << "Matrix Error: " << s << endl;
    cout << flush;
    abort();
    }

void VectorRef::Put0(int i, Real a)
    {
#ifdef MATRIXBOUNDS
    checkindex0(i);
#endif
    if(scale == 0)
	*this = 0.0;
    store[i * stride] = a/scale;
    }

void copyscale(int,double,double*,int,double*,int);
void dcopy(int,double*,int,double*,int);

inline void
VectorRef::assign(const VectorRef &other,Real extrafac)
    {
    extrafac *= other.scale;
    int os = other.stride;
    if(extrafac == 1.0)
	dcopy(length,other.store,os,store,stride);
    else
	copyscale(length,extrafac,other.store,os,store,stride);
    // for (VIter v(*this), o(other); v.test(); v.inc(),o.inc())
// 	v.val() = extrafac * o.val();
    }

VectorRef &
VectorRef::operator = (const VectorRef &other)
    {
    checkcompatibility(other);
    checkassignable();
    if(other.Length() < 1) return *this;
    if(overlaps(*this, other))
	{ *this = Vector (other); return *this; }
    assign(other);
    return *this;
    }

void VectorRef::TreatAsMatrix(MatrixRef& M, int nr, int nc) const
    {
#ifdef MATRIXBOUNDS
    if(stride != 1) error("stride != 1 in TreatAsMatrix");
    if(nr*nc != length) error("nr*nc != length in TreatAsMatrix");
#endif
    M.slink << slink;
    M.nrows = nr;
    M.rowstride = M.ncols = nc;
    M.store = store;
    M.scale = scale;
    M.transpose = 0;
    }

inline void
VectorRef::addin(const VectorRef &other,Real extrafac)
    {
    extrafac *= other.scale;
    if(extrafac != 0.0)
    	{
#if defined(i386) || defined(__x86_64)
	int os = other.stride;
	daxpy_(&length,&extrafac,other.store,&os,store,&stride);
#else
	daxpy(length,extrafac,other.store,other.stride,store,stride); // mine
#endif
    	}
    }

VectorRef &
VectorRef::operator += (const VectorRef &other)
    {
    checkcompatibility(other);
    checkassignable();
    if(other.Length() < 1) return *this;
    if(other.scale == 0.0) return *this;
    if(overlaps(*this,other))
	{ *this += Vector(other); return *this; }
    addin(other);
    return *this;
    }

VectorRef &
VectorRef::operator *= (Real a)
    {
    checkassignable();
    // for (VIter v(*this); v.test(); v.inc())
// 	v.val() *= a;
    void dscal(int n,double a,register double *x,register int incx);
    dscal(length,a,store,stride);
    return *this;
    }

VectorRef &
VectorRef::operator = (Real a)
    {
    checkassignable();
    for (VIter v(*this); v.test(); v.inc())
	v.val() = a;
    return *this;
    }

VectorRef &
VectorRef::operator += (Real a)
    {
    checkassignable();
    for (VIter v(*this); v.test(); v.inc())
	v.val() += a;
    return *this;
    }

void
VectorRef::assign(const MatrixVectorRes & R,int noclear)
    {
    checkassignable();
    if(overlaps(*this,R.a) || overlaps(*this,R.b))
        { 
	if(noclear) *this += (Vector)(R);
	else	    *this  = (Vector)(R);
	return;
        }
    switch(R.op)
        {
        case multiplication:
            mult(R.a,R.b*R.scale,*this,noclear);
            break;
        default:
            _merror("VectorRef: bad R.op in assign MatrixVectorRes");
        }
    }
   
void
VectorRef::assign(const VectorVectorRes & R,int noclear)
    {
    checkassignable();
    if(overlaps(*this,R.a) || overlaps(*this,R.b))
        { 
	if(noclear) *this += (Vector)(R);
	else	    *this  = (Vector)(R);
	return;
        }
	
    switch (R.op)
	{
    case subtraction:
	add(R.a * R.scale, R.b * (-R.scale), *this, noclear);
	break;
    case addition:
	if(R.scale != 1.0) add(R.a * R.scale,	R.b * R.scale,	*this, noclear);
	else		   add(R.a,		R.b,		*this, noclear);
	break;
    default:
	_merror("VectorRef: bad R.op in operator =");
	}
    }

Real
VectorRef::operator * (const VectorRef &other) const
    {
    if (length != other.length)
	_merror("VectorRef *: unequal lengths");
    Real fac = scale * other.scale;
    if(fac == 0.0) return 0.0;
    Real res = 0.0;
    for (VIter v(*this),o(other); v.test(); v.inc(),o.inc())
	res += v.val() * o.val();
    return res * fac;
    }

VectorRef &
VectorRef::operator /= (const VectorRef &other)
    {
    if (length != other.length)
	_merror("VectorRef /=: unequal lengths");
    checkassignable();
    for (VIter v(*this),o(other); v.test(); v.inc(),o.inc())
	v.val() /= other.scale * o.val();
    return *this;
    }

void
add(const VectorRef & A, const VectorRef & B, VectorRef & C,int noclear)
    {
    if (A.length != B.length || A.length != C.length)
	_merror("VectorRef add: unequal lengths");
    C.scale = 1.0;
    if(!noclear) 
	C = A;
    else
	C += A;
    C += B;
    }

Real Norm(const VectorRef &V)
    {
    Real res = 0.0;
    for (VIter v(V); v.test(); v.inc())
	res += v.val()*v.val();
    return fabs(V.scale) * sqrt(res);
    }

Real
VectorRef::sumels() const
    {
    Real res = 0.0;
    for (VIter v(*this); v.test(); v.inc())
	res += v.val();
    return res*scale;
    }

void MatrixRef::Put0(int i,int j, Real a)
    {
    CHECKIND0(i,j);
    if(scale == 0) *this = 0.0;
    store[index0(i,j)] = a/scale;
    }

MatrixRef & 
MatrixRef::operator = (const MatrixRef &other)
    {
    checkassignable();
    checkcompatibility(other);
    if(other.Nrows()*other.Ncols() < 1) return *this;
    if(overlaps(*this,other))
	{ *this = Matrix(other); return *this; }
    RowIter row(*this);
    RowIter orow(other);
    while(row.inc())
	{ orow.inc(); row = orow; }
    return *this;
    }

MatrixRef & 
MatrixRef::operator += (const MatrixRef &other)
    {
    checkassignable();
    checkcompatibility(other);
    if(other.Nrows()*other.Ncols() < 1) return *this;
    if(overlaps(*this,other))
	{ *this += Matrix(other); return *this; }
    RowIter row(*this);
    RowIter orow(other);
    while(row.inc())
	{ orow.inc(); row += orow; }
    return *this;
    }

void
add(const MatrixRef & M1, const MatrixRef & M2, MatrixRef & M3,int noclear)
    {
    if(!noclear)
	M3 = M1;
    else
	M3 += M1;
    M3 += M2;
    }

MatrixRef & 
MatrixRef::operator *= (Real a)
    {
    checkassignable();
    if (nrows > ncols)
	{
	ColumnIter col(*this);
	while(col.inc())
	    col *= a;
	}
    else
	{
	RowIter row(*this);
	while(row.inc())
	    row *= a;
	}
    return *this;
    }

MatrixRef & 
MatrixRef::operator = (Real a)
    {
    checkassignable();
    if (nrows > ncols)
	{
	ColumnIter col(*this);
	while(col.inc())
	    col = 0.0;
	}
    else
	{
	RowIter row(*this);
	while(row.inc())
	    row = 0.0;
	}
    Diagonal() = a;
    return *this;
    }

MatrixRef & 
MatrixRef::operator += (Real a)
    {
    checkassignable();
    Diagonal() += a;
    return *this;
    }

Real MatrixRef::zerofrac() const
    {
    int i, j, nzeros = 0;
    int nr = nrows;
    int nc = ncols;
    for (i = 0; i < nr; i++)
	for (j = 0; j < nc; j++)
	    if (store[i * rowstride + j] == 0.0)
		nzeros++;
    return (Real (nzeros)) /(nr * nc);
    }

void 
mult(const MatrixRef & M, const VectorRef & V, VectorRef & res,int noclear)
    {
    if(!noclear) res = 0.0;
    int i=1;
    ColumnIter mcol(M);
    while(mcol.inc())
	res.addin(mcol,V(i++));
    }

void 
AddOuter(const VectorRef & V1, const VectorRef & V2, MatrixRef & M)
    {
    M.checkassignable();
    if(V1.Length() != M.Nrows() || V2.Length() != M.Ncols())
	_merror("AddOuter: V1.Length != M.Nrows or V2.Length != M.Ncols");
    int i = 1;
    int dobycols = 0;
    dobycols = M.Ncols() < 8 &&  M.Nrows() > M.Ncols() * 2;
    						// Assist caching
    if(dobycols)				// Do by columns
	{
	ColumnIter col(M);
	while(col.inc())
	    col.addin(V1,V2(i++));
	}
    else
	{
	RowIter row(M);
	while(row.inc())
	    row.addin(V2,V1(i++));
	}
    }

#if defined(i386) || defined(__x86_64)
void 
mult(const MatrixRef & M1, const MatrixRef & M2, MatrixRef & M3, int noclear)
    {
#ifdef MATRIXBOUNDS
    if (M1.Ncols() != M2.Nrows())
	_merror("Matrix::mult(M1,M2,M3): Matrices M1, M2 incompatible");
    if (M1.Nrows() != M3.Nrows() || M2.Ncols() != M3.Ncols())
	_merror("Matrix::mult(M1,M2,M3): Matrix M3 incompatible");
#endif

// Use BLAS 3 routine
// Have to reverse the order, since we are really multiplying Ct = Bt*At
    int m = M3.ncols;
    int n = M3.nrows;
    int k = M2.Nrows();
    int lda = M2.rowstride;
    int ldb = M1.rowstride;
    int ldc = M3.rowstride;

    Real beta = noclear ? 1.0 : 0.0;
    Real sca = M1.Scale() * M2.Scale();
    Real *pa = M2.Store();
    Real *pb = M1.Store();
    Real *pc = M3.Store();

    static char pt[] = {'N','T'};
    char transb = pt[M1.DoTranspose()];
    char transa = pt[M2.DoTranspose()];
    dgemm_(&transa,&transb,&m,&n,&k,&sca,pa,&lda,pb,&ldb, &beta, pc, &ldc);
    }

#else

xxxx

void dgemm(const MatrixRef & M1, const MatrixRef & M2, MatrixRef & M3,
		Real alpha, Real beta);
void 
mult(const MatrixRef & M1, const MatrixRef & M2, MatrixRef & M3,int noclear)
    {
    Real beta = noclear ? 1.0 : 0.0;
    Real sca = M1.Scale() * M2.Scale();
    dgemm(M1,M2,M3,sca,beta);		// my own dgemm
    					// It ignores
					// scale factors in M1,M2
    }

/*
void 
mult(const MatrixRef & M1, const MatrixRef & M2, MatrixRef & M3,int noclear)
    {
// M3_ik = M1_ij * M2_jk
// Row, Column, (), Nrows(), Ncols() include tranposes

    if(!noclear) M3 = 0.0;
    int dozeros = 0;	// dozeros: whether to decide by # of zeros
    Real zf1 = 0.0,zf2 = 0.0;
    if(min(M1.Nrows(),M2.Ncols()) > 5)
	{
	dozeros = 1;
	Real minzf = min(zf1 = M1.zerofrac(), zf2 = M2.zerofrac());
	Real maxzf = max(zf1, zf2);
	if(maxzf < 0.7 || minzf/maxzf > 0.7) dozeros = 0;
	}
    if((dozeros && zf1 < zf2) || (!dozeros && M2.Ncols() < M1.Nrows()))
	{					// want M2.Ncols() small
	ColumnIter m2col(M2);
	ColumnIter m3col(M3);
	while(m2col.inc())			// want  zerofrac2 large
	    {
	    m3col.inc();
	    mult(M1,m2col,m3col,1);
	    }
	}
    else
	{
	MatrixRef m2t(M2.t());			// want M1.Nrows() small
	RowIter m1row(M1);
	RowIter m3row(M3);
	while(m1row.inc())				// want  zerofrac1 large
	    {
	    m3row.inc();
	    mult(m2t,m1row,m3row,1);
	    }
	}
    }
*/

#endif

inline Real 
quickran(int & idum)
    {
    const int im = 134456;
    const int ia = 8121;
    const int ic = 28411;
    const Real scale = 1.0 / im;
    idum = (idum*ia+ic)%im;
    return Real(idum) * scale;
    }

VectorRef & 
VectorRef::Randomize()
    {
    static int idum = abs((int) ((long)this));
    quickran(idum);
    for (VIter v(*this); v.test(); v.inc())
	v.val() = quickran(idum);
    return *this;
    }

void
elmult(const VectorRef & V1, const VectorRef & V2,
		      const VectorRef & V3,int noclear)
    {
    V3.checkassignable();
    V1.checkcompatibility(V2);
    V1.checkcompatibility(V3);
    Real scale = V1.scale * V2.scale;
    if(noclear)
	{
	for (VIter v1(V1), v2(V2), v3(V3); v1.test(); v1.inc(),
						v2.inc(),v3.inc())
	    v3.val() += v1.val() * v2.val() * scale;
	}
    else
	{
	for (VIter v1(V1), v2(V2), v3(V3); v1.test(); v1.inc(),
						v2.inc(),v3.inc())
	    v3.val() = v1.val() * v2.val() * scale;
	}
    }

MatrixRef & 
MatrixRef::Randomize()
    {
    if(ncols == 0 || nrows == 0)
	return *this;
    //int nc = Ncols();
    ColumnIter col(*this);
    while(col.inc())
	col.Randomize();
    return *this;
    }

void
MatrixRef::assign(const MatrixMatrixRes & R,int noclear)
    {
    checkassignable();
    if(Nrows() != R.Nrows())
        _merror("MatrixRef assign MatrixmatrixRes: unequal nrows");
    if(Ncols() != R.Ncols())
        _merror("MatrixRef assign MatrixmatrixRes: unequal ncols");
    if (overlaps(*this, R.a) || overlaps(*this, R.b))
	{ 
	if(noclear) *this += (Matrix) (R); 
	else *this = (Matrix) (R); 
	return;
	}

    switch (R.op)
	{
    case multiplication:
	if (R.transpose)
	    mult(R.b.t(), R.a.t() * R.scale, *this, noclear);
	else
	    mult(R.a * R.scale, R.b, *this, noclear);
	break;
    case subtraction:
	if (R.transpose)
	    add(R.a.t() * R.scale, R.b.t() * (-R.scale), *this, noclear);
	else
	    add(R.a * R.scale, R.b * (-R.scale), *this, noclear);
	break;
    case addition:
	if (R.transpose)
	    add(R.a.t() * R.scale, R.b.t() * R.scale, *this, noclear);
	else
	    add(R.a * R.scale, R.b * R.scale, *this, noclear);
	break;
    default:
	_merror("MatrixRef: bad R.op in assign");
	}
    }

ostream & operator << (ostream & s, const VectorRef & V)
    { return s << Vector(V); }
ostream & operator << (ostream & s, const MatrixRef & M)
    { return s << Matrix(M); }
ostream & operator << (ostream & s, const MatrixMatrixRes & R)
    { return s << (Matrix)(R); }
ostream & operator << (ostream & s, const MatrixVectorRes & V)
    { return s << (Vector)(V); }
