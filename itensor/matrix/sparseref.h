// sparseref.h -- Header file for the SparseRef class-- S.R. White 1/00 

#ifndef _sparseref_h
#define _sparseref_h

#include "itensor/matrix/sparse.h"

namespace itensor {

class SparseMatBase : private Array1<SparseVector>
    {

    SparseMatBase(int nr, int nc) 
        : 
        Array1<SparseVector>(nr),
        nrows(nr), 
        ncols(nc), 
        numref(1)
        { }

    SparseVector& 
    Row(int i) 
        { return Array1<SparseVector>::operator[](i); }

    void 
    OwnCopy()
        {
        ArrayDo* a = Array1<SparseVector>::operator->();
        a->OwnCopy();
        }

    friend class SparseRef;

    int nrows; 
    int ncols;
    int numref;
    };

enum ClearFlag { Clear, NoClear};

class SparseRef;

void mult(const SparseRef &, const MatrixRef &,
          MatrixRef &,ClearFlag cf = NoClear);
void mult(const MatrixRef &, const SparseRef &,
          MatrixRef &,ClearFlag cf = NoClear);

class SparseRef
    {
    SparseMatBase * base;
    Real scale;
    int transpose;
    void dodelete()
	{
	if(--base->numref == 0)
	    delete base;
	base = 0;
	}
    SparseVector& Row(int i) const
	{ return base->Row(i); }
    void OwnCopy()
	{
	if(NumRef() > 1)
	    {
	    SparseMatBase* newbase = new SparseMatBase(*base);
	    newbase->OwnCopy();
	    dodelete();
	    base = newbase;
	    base->numref = 1;
	    }
	}
public:
    SparseRef(int nr = 0, int nc = 0) 
	: base(new SparseMatBase(nr,nc)), scale(1.0), transpose(0)
	{}
    SparseRef(const SparseRef & other) 
	: base(other.base), scale(other.scale), transpose(other.transpose)
	{
	base->numref++;
	}
    SparseRef& operator=(const SparseRef &other)
	{
	dodelete();
	base = other.base;
	scale = other.scale;
	transpose = other.transpose;
	base->numref++;
	return *this;
	}
    ~SparseRef() { dodelete(); }
    int Nrows() const { return transpose ? base->ncols : base->nrows; }
    int Ncols() const { return transpose ? base->nrows : base->ncols; }
    int NumRef() const { return base->numref; }
    int RowLen(int row) const { return Row(row).Length(); }
    int Column(int row,int i) const { return Row(row).index(i); }
    Real Element(int row,int i) const { return Row(row).data(i); }
    Real operator() (int i, int j) const
	{
	return scale * (transpose ? Row(j).el(i-1) : Row(i).el(j-1));
	}
    SparseRef operator *(Real a) const
	{
	SparseRef res(*this);
	res.scale *= a;
	return res;
	}
    SparseRef & operator *= (Real a)
	{
	scale *= a;
	return *this;
	}
    SparseRef t() const
	{
	SparseRef res(*this);
	res.transpose = !res.transpose;
	return res;
	}
    inline friend SparseRef operator * (Real, const SparseRef&);

    SparseRef(const MatrixRef&);
    SparseRef& operator+=(const MatrixRef&);
    void PutInMatrix(Matrix&) const;
    friend void mult(const SparseRef &, const MatrixRef &,
		      MatrixRef &,ClearFlag cf);
    friend void mult(const MatrixRef &, const SparseRef &,
		      MatrixRef &,ClearFlag cf);
    void write(std::ostream& s) const;
    void read(std::istream& s);
    void  RemoveElement(int row,int i)
	{
	OwnCopy();
	base->Row(row).data(i) = 0.0;
	}
    void ClearZeroes();
    int memory() const;
    };
inline std::ostream & operator << (std::ostream &s, const SparseRef &a) {return s; }

inline SparseRef operator * (Real a, const SparseRef& sp)
    {
    SparseRef res(sp);
    res.scale *= a;
    return res;   
    }

};

#endif
