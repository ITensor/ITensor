// matrixref.h -- Header file for the MatrixRef class-- S.R. White 8/94 

#ifndef _matrixref_h
#define _matrixref_h

//#define MATRIXBOUNDS		/* Define this for bounds checking  * /

#include "minmax.h"
#include "storelink.h"
#include "tarray1.h"
#include <Accelerate/Accelerate.h>
#include <cassert>

inline void trade(int & i, int & j) { int k = i; i = j; j = k; }

void _merror(const char *);		// error message and exit 

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

#include "matrixref.ih"

#ifdef HEADER_DEFS
ARRAY1CC_DEFS(MatrixRef)
ARRAY1CC_DEFS(VectorRef)

#else //ifndef HEADER_DEFS

#ifdef THIS_IS_MAIN
ARRAY1CC_DEFS(MatrixRef)
ARRAY1CC_DEFS(VectorRef)
#endif

#endif //HEADER_DEFS

#endif
