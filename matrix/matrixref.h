// matrixref.h -- Header file for the MatrixRef class-- S.R. White 8/94 

#ifndef _matrixref_h
#define _matrixref_h

//#define MATRIXBOUNDS		/* Define this for bounds checking  * /

#include "minmax.h"
#include "storelink.h"
#include "tarray1.h"
#include <cassert>
#include <cstdlib>

inline void trade(int & i, int & j) { int k = i; i = j; j = k; }

void _merror(const char *);		// error message and exit 

class MatrixError : public ITError
    {
public:
    typedef ITError
    Parent;

    MatrixError(const std::string& message) 
        : Parent(message)
        { }
    };

class Vector;
class MatrixRef;
class MatrixVectorRes;
class VectorVectorRes;
class VectorRefBare;
class MatrixRefBare;
class VectorRef;
void elmult(const VectorRef &, const VectorRef &, const VectorRef &,int noclear = 0);
void mult(const MatrixRef &, const MatrixRef &, MatrixRef &,int);
void mult(const MatrixRef &, const VectorRef &, VectorRef &,int);
void add(const VectorRef &, const VectorRef &, VectorRef &,int noclear = 0);

class VectorRef		// This class never allocates storage for itself !!
    {
public:
    inline VectorRef SubVector(int l, int u) const;
    inline VectorRef SubVector(int first, int length, int str) const;
    inline VectorRef SubVector0(int l, int u) const;
    inline VectorRef SubVector0(int first, int length, int str) const;

// The following operations have the VectorRef on the left of an =.
// They carry out their actions on the vector referred to.

    VectorRef & operator  = (const VectorRef &);
    VectorRef & operator += (const VectorRef &);
    VectorRef & operator /= (const VectorRef &);	// el by el division
    inline VectorRef & operator -= (const VectorRef &);

    VectorRef & operator = (Real);
    VectorRef & operator *= (Real);
    VectorRef & operator += (Real);
    inline VectorRef & operator /= (Real);

    inline virtual VectorRef & operator = (const VectorVectorRes &);
    inline VectorRef & operator += (const VectorVectorRes &);
    inline VectorRef & operator -= (const VectorVectorRes &);

    inline virtual VectorRef & operator = (const MatrixVectorRes &);
    inline VectorRef & operator += (const MatrixVectorRes &);
    inline VectorRef & operator -= (const MatrixVectorRes &);

// The Ref classes do not allow assignment by element using () or el().
// This is because of difficulties/inefficiencies associated with
// scale and transpose.
// Access to elements is inefficient in general for the Ref classes,
// and should be avoided.
// If you want to assign element by element, for an object with
// scale = 1 and transpose = 0, use the VectorRefBare and
// MatrixRefBare classes, or do it on a Vector or Matrix.

    inline Real operator() (int) const;
    inline Real el(int) const;		// Start from 0
    void Put0(int,Real);
    inline void Put(int,Real);

// Operations that appear to the right of an = return a VectorRef

    inline VectorRef operator * (Real) const;
    inline friend VectorRef operator * (Real, const VectorRef&);
    inline VectorRef operator / (Real fac) const;
    inline VectorRef operator - () const;

// Operations returning a Real;

    Real operator * (const VectorRef&) const;	// dot product
    Real norm() const;
    friend Real Norm(const VectorRef &);	// sqrt(V*V)
    Real sumels() const;			// Sum els of Vector

// Operations on the VectorRefs themselves:

    inline VectorRef & operator<<(const VectorRef &);	// Copy Ref, not vector
    inline VectorRef();
    inline VectorRef(const StoreLink &,Real *,int len,int str=1,Real sca=1.0);
    inline VectorRef(const VectorRef &);

    inline Real* Store() const;			// Allow access to store 
//    inline Real*& AccessStore();	// REALLY allow access to store
    void ShiftStore(int s) 		// potentially dangerous!
	{ store += s; }
    inline int Length() const;
    int& AccessLength()
	{ return length; }		// REALLY allow access to length
    inline int Stride() const;
    inline Real Scale() const;
    inline Real* Last() const;
    inline Real* First() const;

// Length = nrows * ncols
    void TreatAsMatrix(MatrixRef&,int nr, int nc) const;	
// More complicated miscellaneous operations:

    friend void AddOuter(const VectorRef &, const VectorRef &, MatrixRef &);
			// add outer product of vectors into matrix

    VectorRef & Randomize();	
		// Make the elements uniformly random between 0 and 1

    friend void elmult(const VectorRef &, const VectorRef &, const VectorRef &,int noclear);
		// Element by element multiply.
    friend void mult(const MatrixRef &, const MatrixRef &, MatrixRef &,int);
    friend void mult(const MatrixRef &, const VectorRef &, VectorRef &,int);
    friend void add(const VectorRef &, const VectorRef &, VectorRef &,int noclear);

    void assign(const VectorVectorRes &,int);
    void assign(const MatrixVectorRes &,int);

    void write(std::ostream& s) const;

    virtual ~VectorRef () {}

    friend class Vector;
    ARRAY1H_DEFS(VectorRef)
protected:
    Real scale;				// Extra factor applied to vector.
    Real* store;
    int length;
    int stride;
    StoreLink slink;			// Link to original storage, the
					// storage of a Matrix or Vector.
    inline void init();
    inline void copyvars(const VectorRef &V);
    inline void checkcompatibility(const VectorRef &other) const;
    inline void checkassignable() const;
    inline void checkindex(int) const;
    inline void checkindex0(int) const;
    inline void assign(const VectorRef &,Real = 1.0);
    inline void addin(const VectorRef &,Real = 1.0);
    };

Real Norm(const VectorRef &);	// sqrt(V*V)

class VectorRefBare : private VectorRef		// for fast element assignment
    {
public:
    VectorRefBare(const VectorRef &);		// Must have scale == 1.0
    Real  operator () (int) const;
    Real& operator () (int);
    Real  el(int) const;
    Real& el(int);
    };

class Matrix;
class MatrixMatrixRes;
class MatrixRef;
void mult(const MatrixRef &, const MatrixRef &, MatrixRef &,int noclear = 0);
void mult(const MatrixRef &, const VectorRef &, VectorRef &,int noclear = 0);
void add(const MatrixRef &, const MatrixRef &, MatrixRef &,int noclear = 0);

class MatrixRef
    {
public:
// Operations which carry out their actions on the submatrix referred to.
// These have the MatrixRef on the left of an =.

    MatrixRef & operator  = (const MatrixRef &);
    MatrixRef & operator += (const MatrixRef &);
    inline MatrixRef & operator -= (const MatrixRef &);

    MatrixRef & operator = (Real);	// Assign to diagonal of submatrix.
    MatrixRef & operator *= (Real);
    MatrixRef & operator += (Real);	// SubMatrix plus a diagonal Matrix
    inline MatrixRef & operator -= (Real a);

    inline virtual MatrixRef & operator  = (const MatrixMatrixRes &);
    inline MatrixRef & operator += (const MatrixMatrixRes &);
    inline MatrixRef & operator -= (const MatrixMatrixRes &);

// See the note about element access in VectorRef
// Number from start of submatrix.
    inline Real operator() (int, int) const;	// start from 1
    inline Real el(int, int) const;		// Start from 0
    void Put0(int,int,Real a);
    void Put(int,int,Real a);

// Make a new Ref which refers to a part of this one:
    inline VectorRef Diagonal() const;
    inline VectorRef Row(int) const;
    inline VectorRef Column(int) const;
    inline MatrixRef SubMatrix
    	(int firstrow, int lastrow, int firstcol, int lastcol) const;
    inline MatrixRef SubMatrix (int lastrow,int lastcol) const;
    inline MatrixRef Columns(int firstcol, int lastcol) const;
    inline MatrixRef Rows(int firstrow, int lastrow) const;

    inline VectorRef Column0(int) const;
    inline VectorRef Row0(int) const;
    inline MatrixRef SubMatrix0(int, int, int, int) const;
    inline MatrixRef SubMatrix0(int, int) const;
    inline MatrixRef Columns0(int, int) const;
    inline MatrixRef Rows0(int, int) const;

    inline VectorRef TreatAsVector() const;	// Length = nrows * ncols

// Operations that appear to the right of an = return a MatrixRef
    inline MatrixRef operator *(Real) const;
    inline friend MatrixRef operator * (Real, const MatrixRef&);
    inline MatrixRef operator / (Real) const;
    inline MatrixRef operator - () const;	// Negate
    inline MatrixRef t() const;
    inline MatrixRef t(int dotrans) const;	// dotrans = 0, no transpose
    inline friend MatrixRef Transpose(const MatrixRef &);

// Operations on the MatrixRefs themselves
// copy the Ref, not the Matrix
    inline MatrixRef& operator<<(const MatrixRef &);	
    inline MatrixRef(const MatrixRef &);
    inline MatrixRef();

    inline Real* Store() const;		 	// Allow access to store 
    inline int Nrows() const;
    inline int Ncols() const;
    inline int RowStride() const;
    inline Real Scale() const;
    inline int DoTranspose() const;
    inline Real* Last() const;
    inline Real* First() const;
    inline int memory() const;         // return memory used in bytes
    inline int NumRef() const;

// More complicated miscellaneous operations:

    MatrixRef & Randomize();	
		// Make the elements uniformly random between 0 and 1
    inline friend Real Trace(const MatrixRef&);
    friend Matrix Inverse(const MatrixRef &);
    friend void AddOuter(const VectorRef &, const VectorRef &, MatrixRef &);
		// add outer product of vectors into matrix
    Real zerofrac() const;		// Fraction of els which == 0.0

    void write(std::ostream& s) const;

    friend void copybare(MatrixRef &, const MatrixRef &);
    void assign(const MatrixMatrixRes &,int);	// used in operator = and +=
    friend void mult(const MatrixRef &, const MatrixRef &, MatrixRef &,int noclear);
    friend void mult(const MatrixRef &, const VectorRef &, VectorRef &,int noclear);
    friend void add(const MatrixRef &, const MatrixRef &, MatrixRef &,int noclear);

    virtual ~MatrixRef () {}

    friend class Matrix;
    friend class RowIter;
    friend class ColumnIter;
    friend class VectorRef;
    friend class VectorRefNoLink;
    friend class MatrixRefNoLink;
    ARRAY1H_DEFS(MatrixRef)
protected:
    Real *store;
    Real scale;
    StoreLink slink;
    int nrows;
    int ncols;
    int rowstride;		// Spacing in elements between rows
    int transpose;

    inline void checkcompatibility(const MatrixRef &) const;
    inline void checkassignable() const;
    inline void init();
    inline void copyvars(const MatrixRef &);
    inline int index(int, int) const;
    inline int index0(int, int) const;
    inline void checkindex(int,int) const;
    inline void checkindex0(int,int) const;
    };

class MatrixRefBare : private MatrixRef	// for fast element assignment
    {
public:
    inline MatrixRefBare(const MatrixRef &);
    inline Real  operator () (int, int) const;
    inline Real& operator () (int, int);
    inline Real  el(int, int) const;
    inline Real& el(int, int);
    };

enum binop {multiplication, addition, subtraction};

class VectorVectorRes
    {
public:
    inline VectorVectorRes
    	(const VectorRef &, const VectorRef &, binop, Real sca = 1.0);
    inline VectorVectorRes operator *(Real) const;
    inline VectorVectorRes operator /(Real) const;
    inline VectorVectorRes operator -() const;
    inline int Length() const;

    inline VectorVectorRes(const VectorVectorRes &);
    inline operator Vector() const;
public:
    VectorRef a, b;
    binop op;
    Real scale;
    };

class MatrixMatrixRes
    {
public:
    inline MatrixMatrixRes(const MatrixRef &, const MatrixRef &, binop,
			Real sca = 1.0, int tran = 0);
    inline MatrixMatrixRes t() const;
    inline MatrixMatrixRes operator *(Real) const;
    inline MatrixMatrixRes operator /(Real) const;
    inline MatrixMatrixRes operator -() const;

    inline MatrixMatrixRes(const MatrixMatrixRes &);
    inline int Nrows() const;
    inline int Ncols() const;
    inline MatrixMatrixRes();
    inline operator Matrix() const;
public:
    MatrixRef a, b;
    binop op;
    Real scale;
    int transpose;
    };

class MatrixVectorRes
    {
public:
    inline MatrixVectorRes(const MatrixRef &, const VectorRef &, binop,
			    Real sca = 1.0);
    inline MatrixVectorRes operator *(Real) const;
    inline MatrixVectorRes operator /(Real) const;
    inline MatrixVectorRes operator -() const;

    inline MatrixVectorRes(const MatrixVectorRes &);
    inline int Length() const;
    inline MatrixVectorRes();
    inline operator Vector() const;
//     inline friend Vector& vector(const MatrixVectorRes &);
public:
    MatrixRef a;
    VectorRef b;
    binop op;
    Real scale;
    };

// These have to be after the declarations of VectorVectorRes, etc.


class RowIter : public VectorRef
    {
public:
// Usage:
// RowIter row(M);
// while(row.inc()) { ... }
    inline RowIter(const MatrixRef & M);
    inline int inc();

    inline VectorRef & operator = (const RowIter &);
    inline VectorRef & operator = (const VectorRef &);
    inline VectorRef & operator = (const VectorVectorRes &);
    inline VectorRef & operator = (const MatrixVectorRes &);
    inline VectorRef & operator = (Real);
	
public:
    int rowinc;
    Real* finish;
    };

class ColumnIter : public VectorRef
    {
public:
    inline ColumnIter(const MatrixRef & M);
    inline int inc();

    inline VectorRef & operator = (const ColumnIter &);
    inline VectorRef & operator = (const VectorRef &);
    inline VectorRef & operator = (const VectorVectorRes &);
    inline VectorRef & operator = (const MatrixVectorRes &);
    inline VectorRef & operator = (Real a);

public:
    int colinc;
    Real* finish;
    };

class VIter
    {
public:
// Usage:
// for(VIter viter(V); viter.test(); viter.inc())
//	viter.val() = ...

    inline VIter(const VectorRef &);
    inline void inc();
    inline int test() const;
    inline Real & val();	// Note: this does not give/use scale: where
    				// needed, scale must be put in by hand:
				// ... = viter.val() * V.scale ...
	
public:
    int stride;
    Real *store;
    Real *finish;
    };

class MatrixRefNoLink : public MatrixRef	// for fast Ref operations
    {
public:
    MatrixRefNoLink()  {}
    void operator<<(const MatrixRef &M)
	{
	store = M.store; nrows = M.nrows; ncols = M.ncols;
	rowstride = M.rowstride; transpose = M.transpose; 
	scale = M.scale;
	}
    void SetScale(Real a)
	{ scale = a; }
    void ApplyTrans()
	{ transpose = !transpose; }
    };

class VectorRefNoLink : public VectorRef	// for fast Ref operations
    {
public:
    VectorRefNoLink()  {}
    void operator<<(const VectorRef &V)
	{
	store=V.Store(); length=V.Length(); stride=V.Stride(); scale=V.Scale();
	}
    void operator<<(const MatrixRef &M)
	{
	if(M.rowstride != M.ncols)
	    _merror("bad call to VectorRefNoLink<<MatrixRef");
	store=M.store; length=M.nrows*M.ncols; stride=1; scale=M.scale;
	}
    void SetScale(Real a)
	{ scale = a; }
    VectorRef & operator = (Real a)
	{ return VectorRef::operator=(a); }
    VectorRef & operator = (const VectorRef &other)
	{ return VectorRef::operator=(other); }
    VectorRef & operator = (const VectorVectorRes &other)
	{ return VectorRef::operator=(other); }
    VectorRef & operator = (const MatrixVectorRes &other)
	{ return VectorRef::operator=(other); }
    };

std::ostream & operator << (std::ostream &s, const MatrixRef &a);
std::ostream & operator << (std::ostream &s, const VectorRef &a);
std::ostream & operator << (std::ostream &s, const MatrixMatrixRes &a);
std::ostream & operator << (std::ostream &s, const MatrixVectorRes &a);

// Inlines for matrixref.h
#ifdef MATRIXBOUNDS
#define	CHECKINDEX(i)	checkindex(i)
#define	CHECKINDEX0(i)	checkindex0(i)
#define	CHECKIND(i,j)	checkindex(i,j)
#define	CHECKIND0(i,j)	checkindex0(i,j)
#else
#define	CHECKINDEX(i)
#define	CHECKINDEX0(i)
#define	CHECKIND(i,j)
#define	CHECKIND0(i,j)
#endif

inline Real* VectorRef::Store() const { return store; }			

// inline Real*& VectorRef::AccessStore() { return store; }			

inline int VectorRef::Length() const { return length; }

inline int VectorRef::Stride() const { return stride; }

inline Real VectorRef::Scale() const { return scale; }

inline Real* VectorRef::Last() const { return store+(length-1)*stride; }

inline Real* VectorRef::First() const { return store; }

inline VectorRef::VectorRef() { init();}

inline VectorRef::VectorRef(const StoreLink & sl,Real *st, 
    int len, int str, Real sca) : slink(sl) 
    { store=st; length=len; stride=str; scale=sca; }

inline VectorRef::VectorRef(const VectorRef & V) 
    : slink(V.slink) { copyvars(V); }

inline VectorRef & VectorRef::operator<<(const VectorRef & V)		
    { slink << V.slink; copyvars(V); return *this; }

inline Real VectorRef::operator() (int i) const	
    { CHECKINDEX(i); return scale * store[(i - 1) * stride]; }

inline Real VectorRef::el(int i) const		
    { CHECKINDEX0(i); return scale * store[i * stride]; }

inline void VectorRef::Put(int i,Real a) { Put0(i-1,a); }

inline VectorRef & VectorRef::operator-= (const VectorRef & other)
    { *this += -other; return *this; }

inline VectorRef & VectorRef::operator= (const VectorVectorRes &R)
    { assign(R,0); return *this; }

inline VectorRef & VectorRef::operator+= (const VectorVectorRes &R)
    { assign(R,1); return *this; }

inline VectorRef & VectorRef::operator= (const MatrixVectorRes &R)
    { assign(R,0); return *this; }

inline VectorRef & VectorRef::operator+= (const MatrixVectorRes &R)
    { assign(R,1); return *this; }

inline VectorRef & VectorRef::operator/= (Real fac)
    { *this *= 1.0/fac; return *this; }

inline VectorRef VectorRef::operator* (Real fac) const
    { return VectorRef(slink,store,length,stride,scale*fac); } 

inline VectorRef operator* (Real fac, const VectorRef& V)
    { return V * fac; }

inline VectorRef VectorRef::operator/ (Real fac) const
    { return VectorRef(slink,store,length,stride,scale/fac); } 

inline VectorRef VectorRef::operator- () const
    { return VectorRef(slink,store,length,stride,-scale); } 

inline VectorRef VectorRef::SubVector(int l, int u) const
    { 
    if(u > length || l < 1 || u < l)
	_merror("bad call to SubVector");
    return VectorRef(slink,store+(l-1)*stride, u-l+1, stride, scale); 
    }

inline VectorRef VectorRef::SubVector0(int l, int u) const
    {
    return SubVector(l+1,u+1);
    }

inline VectorRef VectorRef::SubVector(int first, int len, int str) const
    { 
    if(length <= 0 || first < 1 || first + (len-1) * str > length)
	_merror("bad call to SubVector");
    return VectorRef(slink,store+(first-1)*stride, len, stride * str, scale); 
    }

inline VectorRef VectorRef::SubVector0(int first, int len, int str) const
    {
    return SubVector(first+1,len,str);
    }

inline void VectorRef::init()
    { store = 0; length = 0; stride = 0; scale = 1.0; }

inline void VectorRef::copyvars(const VectorRef &V) 
    { store=V.store; length=V.length; stride=V.stride; scale=V.scale; }

inline void VectorRef::checkcompatibility(const VectorRef &other) const
    { if (length != other.length)_merror("VectorRef: unequal lengths");}

inline void VectorRef::checkassignable() const
    { if (scale != 1.0) _merror("VectorRef assignment: scale != 1.0"); }

inline void VectorRef::checkindex(int i) const
    {
    if (i < 1 || i > length)
	{
    std::cerr << "index=(" << i << ")\n";
	_merror("VectorRef: index out of bounds");
	}
    }

inline void VectorRef::checkindex0(int i) const
    { checkindex(i+1); }

inline VectorRef & VectorRef::operator -= (const VectorVectorRes &R)
    { *this += -R; return *this; }

inline VectorRef & VectorRef::operator -= (const MatrixVectorRes &R)
    { *this += -R; return *this; }

inline VectorVectorRes 
operator + (const VectorRef & V1, const VectorRef &V2)
    { return VectorVectorRes(V1,V2,addition); }

inline VectorVectorRes 
operator - (const VectorRef & V1, const VectorRef &V2)
    { return VectorVectorRes(V1,-V2,addition); }

inline MatrixVectorRes 
operator * (const MatrixRef & M, const VectorRef &V)
    { return MatrixVectorRes(M,V,multiplication); }

inline MatrixVectorRes 
operator * (const VectorRef &V, const MatrixRef & M)
    { return MatrixVectorRes(M.t(),V,multiplication); }

inline MatrixMatrixRes 
operator + (const MatrixRef & M1, const MatrixRef &M2)
    { return MatrixMatrixRes(M1,M2,addition); }

inline MatrixMatrixRes 
operator - (const MatrixRef & M1, const MatrixRef &M2)
    { return MatrixMatrixRes(M1,-M2,addition); }

inline MatrixMatrixRes 
operator * (const MatrixRef & M1, const MatrixRef &M2)
    { return MatrixMatrixRes(M1,M2,multiplication); }

inline VectorVectorRes 
operator * (Real a, const VectorVectorRes & A)
    { return A * a; }

inline MatrixMatrixRes 
Transpose(const MatrixMatrixRes & A)
    { return A.t(); }

inline MatrixMatrixRes 
operator * (Real a, const MatrixMatrixRes & A)
    { return A * a; }

inline MatrixVectorRes 
operator * (Real a, const MatrixVectorRes & A)
    { return A * a; }

// Two VectorRefs are free from overlap if the are in separate
// locations in memory, or if the have identical strides but
// are offset by a non-multiple of the stride

inline int 
overlaps(const VectorRef & A, const VectorRef & B)
    { return !( A.Last() < B.First() || B.Last() < A.First() 
    || (A.Stride() == B.Stride() && (A.Last()-B.First())%A.Stride() != 0 ) ); }

inline int 
overlaps(const MatrixRef & A, const MatrixRef & B)
    { return !( A.Last() < B.First() || B.Last() < A.First() ); }

inline int 
overlaps(const MatrixRef & A, const VectorRef & B)
    { return !( A.Last() < B.First() || B.Last() < A.First() ); }

inline int 
overlaps(const VectorRef & A, const MatrixRef & B)
    { return overlaps(B,A); }

inline VectorRefBare::VectorRefBare(const VectorRef & V) : VectorRef(V)
    { if(scale != 1.0) _merror("scale != 1.0 in VectorRefBare"); }

inline Real VectorRefBare::operator () (int i) const
    { CHECKINDEX(i); return store[(i-1) * stride]; }

inline Real & VectorRefBare::operator () (int i) 
    { CHECKINDEX(i); return store[(i-1) * stride]; }

inline Real   VectorRefBare::el(int i) const
    { CHECKINDEX0(i); return store[i * stride]; }

inline Real & VectorRefBare::el(int i) 
    { CHECKINDEX0(i); return store[i * stride]; }

inline Real* MatrixRef::Store() const
    { return store; }			// Allow access to store 

inline int MatrixRef::Nrows() const
    { return (transpose ? ncols : nrows); }

inline int MatrixRef::Ncols() const
    { return (transpose ? nrows : ncols); }

inline int MatrixRef::RowStride() const
    { return rowstride; }

inline Real MatrixRef::Scale() const
    { return scale; }

inline int MatrixRef::DoTranspose() const
    { return transpose; }

inline Real* MatrixRef::Last() const
    { return store+(nrows-1)*rowstride+ncols-1; }

inline Real* MatrixRef::First() const
    { return store; }

inline MatrixRef::MatrixRef () 
    { init(); }

inline MatrixRef::MatrixRef(const MatrixRef &M) : slink(M.slink)
	{ copyvars(M); }

inline MatrixRef& MatrixRef::operator<<(const MatrixRef & M)
    { slink << M.slink; copyvars(M); return *this; }

inline Real MatrixRef::operator() (int i, int j) const 
    { CHECKIND(i,j); return scale * store[index(i,j)]; }

inline Real MatrixRef::el(int i, int j) const
    { CHECKIND0(i,j); return scale * store[index0(i,j)]; }

inline void MatrixRef::Put(int i,int j,Real a)
    { Put0(i-1,j-1,a); }

inline VectorRef MatrixRef::Row(int i) const
    {
    if(i>Nrows() || i<1) _merror("VectorRef::Row i > Nrows() || i < 1");
    return VectorRef(slink, store + (i - 1) * (transpose ? 1 : rowstride),
	     transpose ? nrows : ncols, transpose ? rowstride : 1,scale);
    }

inline VectorRef MatrixRef::Column(int i) const
    {
    if(i>Ncols() || i<1)_merror("VectorRef::Column i > Ncols() || i < 1");
    return VectorRef(slink, store + (i - 1) * (transpose ? rowstride : 1),
	 transpose ? ncols : nrows, transpose ? 1 : rowstride,scale);
    }

inline VectorRef MatrixRef::Diagonal() const
    { return VectorRef(slink,store,min(nrows,ncols),rowstride+1,scale); }

inline MatrixRef MatrixRef::SubMatrix(int il, int iu, int jl, int ju) const
    { 
    if(iu > Nrows() || il < 1 || iu < il) _merror("bad call to SubMatrix");
    if(ju > Ncols() || jl < 1 || ju < jl) _merror("bad call to SubMatrix");
    if(transpose) { trade(il,jl); trade(iu,ju); }
    MatrixRef res(*this);
    res.nrows = iu - il + 1;
    res.ncols = ju - jl + 1;
    res.store += (il - 1) * rowstride + jl - 1;
    return res;
    }

inline MatrixRef MatrixRef::SubMatrix(int iu, int ju) const
    { return SubMatrix(1, iu, 1, ju); }

inline MatrixRef MatrixRef::Columns(int jl, int ju) const
    { return SubMatrix(1,Nrows(),jl,ju); }

inline MatrixRef MatrixRef::Rows(int il, int iu) const
    { return SubMatrix(il,iu,1,Ncols()); }

inline VectorRef MatrixRef::Column0(int i) const
    { return Column(i+1); }

inline VectorRef MatrixRef::Row0(int i) const
    { return Row(i+1); }

inline MatrixRef MatrixRef::SubMatrix0(int il, int iu, int jl, int ju) const
    { return SubMatrix(il+1,iu+1,jl+1,ju+1); }

inline MatrixRef MatrixRef::SubMatrix0(int iu,int ju) const
    { return SubMatrix(1,iu+1,1,ju+1); }

inline MatrixRef MatrixRef::Columns0(int jl, int ju) const
    { return Columns(jl+1,ju+1); }

inline MatrixRef MatrixRef::Rows0(int il, int iu) const
    { return Rows(il+1,iu+1); }

inline VectorRef MatrixRef::TreatAsVector() const
    {
    if(rowstride != ncols)
	_merror("bad call to TreatAsVector");
    return VectorRef(slink, store, nrows * ncols, 1, scale);
    }

inline MatrixRef & MatrixRef::operator -= (const MatrixRef & other)
    { *this += -other; return *this; }

inline MatrixRef & MatrixRef::operator  = (const MatrixMatrixRes &R)
    { assign(R,0); return *this; }

inline MatrixRef & MatrixRef::operator += (const MatrixMatrixRes &R)
    { assign(R,1); return *this; }

inline MatrixRef & MatrixRef::operator -= (const MatrixMatrixRes &R)
    { MatrixMatrixRes rr(R); rr.scale *= -1.0; assign(rr,1); return *this; }

inline MatrixRef & MatrixRef::operator -= (Real a)
    { *this += -a; return *this; }

inline MatrixRef MatrixRef::t() const
    { MatrixRef res(*this); res.transpose = !res.transpose; return res; }

inline MatrixRef MatrixRef::t(int dotrans) const
    {
    MatrixRef res(*this);
    if (dotrans)
	res.transpose = !res.transpose;
    return res;
    }

inline MatrixRef Transpose(const MatrixRef & M) 
    { return M.t(); }

inline MatrixRef MatrixRef::operator *(Real fac) const
    { MatrixRef res(*this); res.scale *= fac; return res; }

inline MatrixRef operator* (Real fac, const MatrixRef& V)
    { return V * fac; }

inline MatrixRef MatrixRef::operator/ (Real fac) const
    { return *this * (1.0/fac); }

inline MatrixRef MatrixRef::operator- () const	// Negate 
    { return *this * (-1.0); }

inline Real Trace(const MatrixRef& M)
    { return M.Diagonal().sumels(); }

inline void MatrixRef::checkcompatibility(const MatrixRef &other) const
    {
    if(Nrows() != other.Nrows())
	_merror("MatrixRef compatibility: unequal nrows");
    if(Ncols() != other.Ncols())
	_merror("MatrixRef compatibility: unequal ncols");
    }

inline void MatrixRef::checkassignable() const
    {
    if (transpose != 0)
	_merror("MatrixRef assignment: transpose != 0");
    if (scale != 1.0)
	_merror("MatrixRef assignment: scale != 1.0");
    }

inline void MatrixRef::init() 
    { store = 0; scale = 1.0; nrows = ncols = rowstride = transpose = 0; }

inline void MatrixRef::copyvars(const MatrixRef &M) 
    {  
    store = M.store; nrows = M.nrows; ncols = M.ncols;
    rowstride = M.rowstride; transpose = M.transpose; scale = M.scale;
    }

inline int MatrixRef::index0(int i, int j) const
    { return transpose ? j * rowstride + i : i * rowstride + j; }

inline int MatrixRef::index(int i, int j) const
    { return index0(i-1,j-1); }

inline void MatrixRef::checkindex(int i,int j) const
    {
    if (i < 1 || i > Nrows() || j < 1 || j > Ncols())
	{
    std::cerr << "index=(" << i << "," << j << ") ";
	_merror("MatrixRef: index out of bounds");
	}
    }

inline void MatrixRef::checkindex0(int i,int j) const
    { checkindex(i+1,j+1); }

inline int MatrixRef::memory() const
    { return sizeof(MatrixRef) + slink.memory()/slink.NumRef(); }

inline int MatrixRef::NumRef() const
    { return slink.NumRef(); }

inline Real  MatrixRefBare::operator () (int i, int j) const
    { CHECKIND(i,j); return store[(i-1) * rowstride + j - 1]; }

inline Real& MatrixRefBare::operator () (int i, int j) 
    { CHECKIND(i,j); return store[(i-1) * rowstride + j - 1]; }

inline Real  MatrixRefBare::el (int i, int j) const
    { CHECKIND0(i,j); return store[i * rowstride + j]; }

inline Real& MatrixRefBare::el( int i, int j) 
    { CHECKIND0(i,j); return store[i * rowstride + j]; }

inline MatrixRefBare::MatrixRefBare(const MatrixRef & V)
    	: MatrixRef(V)
    {
    if(scale != 1.0) _merror("scale != 1.0 in MatrixRefBare");
    if(transpose) _merror("transpose != 0 in MatrixRefBare");
    }

inline int VectorVectorRes::Length() const
    { return a.Length(); }

inline VectorVectorRes::VectorVectorRes(const VectorVectorRes & A) 
    : a(A.a), b(A.b) { op = A.op; scale = A.scale; }

inline VectorVectorRes::VectorVectorRes
	(const VectorRef & ax, const VectorRef & bx, binop opx,
	    Real sca) : a(ax), b(bx)
    {
    if (ax.Length() != bx.Length())
	_merror("ax.Length() != bx.Length() in VectorVectorRes");
    if (ax.Length() == 0)
	_merror("ax.Length() == 0 in VectorVectorRes");
    op = opx;
    scale = sca;
    }

inline VectorVectorRes VectorVectorRes::operator *(Real aa) const
    { VectorVectorRes res(*this); res.scale *= aa; return res; }

inline VectorVectorRes VectorVectorRes::operator / (Real fac) const
    { return *this * (1.0/fac); }

inline VectorVectorRes VectorVectorRes::operator - () const
    { return (*this) * -1.0; }

inline int MatrixMatrixRes::Nrows() const
    { return transpose ? b.Ncols() : a.Nrows(); }

inline int MatrixMatrixRes::Ncols() const
    { return transpose ? a.Nrows() : b.Ncols(); }

inline MatrixMatrixRes::MatrixMatrixRes() { }

inline MatrixMatrixRes::MatrixMatrixRes(const MatrixMatrixRes & A) 
	: a(A.a), b(A.b)
    { op = A.op; scale = A.scale; transpose = A.transpose; }

inline MatrixMatrixRes::MatrixMatrixRes(const MatrixRef &ax, 
	const MatrixRef &bx, binop opx, Real sca, int tran) 
	: a(ax), b(bx)
    {
    op = opx;
    scale = sca;
    transpose = tran;
    }

inline MatrixMatrixRes MatrixMatrixRes::t() const
    {
    MatrixMatrixRes res(*this);
    res.transpose = !res.transpose;
    return res; 
    }

inline MatrixMatrixRes MatrixMatrixRes::operator *(Real fac) const
    { MatrixMatrixRes res(*this); res.scale *= fac; return res; }

inline MatrixMatrixRes MatrixMatrixRes::operator / (Real fac) const
    { return *this * (1.0/fac); }

inline MatrixMatrixRes MatrixMatrixRes::operator - () const
    { return (*this) * -1.0; }

inline int MatrixVectorRes::Length() const
    { return a.Nrows(); }

inline MatrixVectorRes::MatrixVectorRes() { }

inline MatrixVectorRes::MatrixVectorRes(const MatrixVectorRes & A) 
    : a(A.a), b(A.b) { op = A.op; scale = A.scale; }

inline MatrixVectorRes::MatrixVectorRes(const MatrixRef &ax, 
	const VectorRef & bx, binop opx, Real sca) 
    : a(ax), b(bx) 
    {
    if(ax.Ncols() != bx.Length())
	_merror("ax.Ncols() != bx.Length()) in MatrixVectorRes");
    op = opx;
    scale = sca;
    }

inline MatrixVectorRes MatrixVectorRes::operator *(Real aa) const
    { MatrixVectorRes res(*this); res.scale *= aa; return res; }

inline MatrixVectorRes MatrixVectorRes::operator / (Real fac) const
    { return *this * (1.0/fac); }

inline MatrixVectorRes MatrixVectorRes::operator - () const
    { return (*this) * -1.0; }

inline RowIter::RowIter(const MatrixRef & M) : VectorRef(M.Row(1))
    {
    rowinc = M.index(2, 1);
    finish = store + M.Nrows() * rowinc;
    store -= rowinc;
    }

inline int RowIter::inc()
    { store += rowinc; return store < finish; }

inline VectorRef & RowIter::operator = (const VectorRef &other)
    { return VectorRef::operator=(other); }

inline VectorRef & RowIter::operator = (const RowIter &other)
    { return VectorRef::operator=(other); }

inline VectorRef & RowIter::operator = (const VectorVectorRes &other)
    { return VectorRef::operator=(other); }

inline VectorRef & RowIter::operator = (const MatrixVectorRes &other)
    { return VectorRef::operator=(other); }

inline VectorRef & RowIter::operator = (Real a)
    { return VectorRef::operator=(a); }

inline ColumnIter::ColumnIter(const MatrixRef & M) 
	: VectorRef(M.Column(1))
    {
    colinc = M.index(1, 2);
    finish = store + M.Ncols() * colinc;
    store -= colinc;
    }

inline int ColumnIter::inc()
    { store += colinc; return store < finish; }

inline VectorRef & ColumnIter::operator = (const VectorRef &other)
    { return VectorRef::operator=(other); }

inline VectorRef & ColumnIter::operator = (const ColumnIter &other)
    { return VectorRef::operator=(other); }

inline VectorRef & ColumnIter::operator = (const VectorVectorRes &other)
    { return VectorRef::operator=(other); }

inline VectorRef & ColumnIter::operator = (const MatrixVectorRes &other)
    { return VectorRef::operator=(other); }

inline VectorRef & ColumnIter::operator = (Real a)
    { return VectorRef::operator=(a); }

inline VIter::VIter(const VectorRef & V)
    {
    stride = V.Stride(); store = V.Store();
    finish = store + V.Length() * stride;
    }

inline void VIter::inc()
    { store += stride; }

inline int VIter::test() const
    { return store < finish; }

inline Real & VIter::val()
    { return *store; }


ARRAY1CC_DEFS(MatrixRef)
ARRAY1CC_DEFS(VectorRef)

#endif
