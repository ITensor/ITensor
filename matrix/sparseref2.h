// sparseref.h -- Header file for the SparseRef class-- S.R. White 2/00 

#ifndef _sparseref_h
#define _sparseref_h

//#include "sparse.h"
#include <map>

namespace itensor {

typedef map<int,Real>::const_iterator CRowIter;

class SparseVec : public map<int,Real>
    {
    SparseVec() {}
    int Storage() const { return size() * (2*sizeof(int) + sizeof(Real)); }
					// estimate
    Real operator()(int i) const
	{
	iterator p = find(i);
	if(p != end())
	    return p->second;
	return 0.0;
	}
    void write(ostream& s) const;
    void read(istream& s);

    };

class SparseMatBase : private Array1<SparseVec>
    {
    friend class SparseRef;
    int nrows; 
    int ncols;
    int numref;
    SparseMatBase(int nr, int nc) 
	: nrows(nr), ncols(nc), numref(1), Array1<SparseVec>(nr) { }
    SparseVec& Row(int i) 
	{ return Array1<SparseVec>::operator[](i); }
    void OwnCopy()
	{
	ArrayDo* a = Array1<SparseVec>::operator->();
	a->OwnCopy();
	}
    };

enum ClearFlag { Clear, NoClear};
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
    SparseVec& Row(int i) const
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
    SparseRef::~SparseRef() { dodelete(); }
    int Nrows() const { return transpose ? base->ncols : base->nrows; }
    int Ncols() const { return transpose ? base->nrows : base->ncols; }
    int NumRef() const { return base->numref; }
//    int RowLen(int row) const { return Row(row).Length(); }
//    int Column(int row,int i) const { return Row(row).index(i); }
//    Real Element(int row,int i) const { return Row(row).data(i); }
    Real operator() (int i, int j) const
	{
	return scale * (transpose ? Row(j)(i) : Row(i)(j));
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
		      MatrixRef &,ClearFlag cf = NoClear);
    friend void mult(const MatrixRef &, const SparseRef &,
		      MatrixRef &,ClearFlag cf = NoClear);
    void write(ostream& s) const;
    void read(istream& s);
/*
    void  RemoveElement(int row,int i)
	{
	OwnCopy();
	base->Row(row).data(i) = 0.0;
	}
    void ClearZeroes();
*/
    int memory() const;
    };
inline ostream & operator << (ostream &s, const SparseRef &a) {return s; }

inline SparseRef operator * (Real a, const SparseRef& sp)
    {
    SparseRef res(sp);
    res.scale *= a;
    return res;   
    }


};

#endif
