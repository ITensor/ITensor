#include "indexset.h"
using namespace std;
using boost::format;
using boost::array;

IndexSet::
IndexSet()
    :
    rn_(0),
    r_(0),
    ur_(0)
    { }

IndexSet::
IndexSet(const Index& i1)
    :
    rn_((i1.m() == 1 ? 0 : 1)),
    r_(1),
    ur_(i1.uniqueReal())
    { 
    index_[1] = i1;
    }

IndexSet::
IndexSet(const Index& i1, const Index& i2)
    :
    r_(2),
    ur_(i1.uniqueReal() + i2.uniqueReal())
    { 
	if(i1.m()==1) 
	    {
	    index_[1] = i2; 
        index_[2] = i1; 
	    rn_ = (i2.m() == 1 ? 0 : 1);
	    }
	else 
	    { 
	    index_[1] = i1; 
        index_[2] = i2; 
	    rn_ = (i2.m() == 1 ? 1 : 2); 
	    }
    }

IndexSet::
IndexSet(Index i1, Index i2, Index i3,
         Index i4, Index i5, Index i6,
         Index i7, Index i8)
    :
    r_(3)
    { 
	array<Index,NMAX> ii = {{ i1, i2, i3, i4, i5, i6, i7, i8 }};
	while(ii[r_] != Index::Null()) ++r_;
    int alloc_size;
    sortIndices(ii,r_,rn_,alloc_size,index_,0);
    }

IndexSet::
IndexSet(const IndexSet& other, const Permutation& P)
    :
    rn_(other.rn_),
    r_(other.r_),
    ur_(other.ur_)
    {
    for(int j = 1; j <= r_; ++j)
        index_[P.dest(j)] = other.index_[j];
    }

Index IndexSet::
findtype(IndexType t) const
	{
    for(int j = 1; j <= rn_; ++j)
    if(index_[j].type() == t) return index_[j];
    Error("IndexSet::findtype failed."); return Index();
	}

bool IndexSet::
findtype(IndexType t, Index& I) const
	{
    for(int j = 1; j <= r_; ++j)
    if(index_[j].type() == t)
        {
        I = index_[j];
        return true;
        }
    return false;
	}

int IndexSet::
findindex(const Index& I) const
    {
    if(I.m() == 1) return findindex1(I);
    else           return findindexn(I);
    return 0;
    }

int IndexSet::
findindexn(const Index& I) const
	{
    for(int j = 1; j <= rn_; ++j)
    if(index_[j] == I) return j;
    return 0;
	}

int IndexSet::
findindex1(const Index& I) const
	{
    for(int j = rn_+1; j <= r_; ++j)
    if(index_[j] == I) return j;
    return 0;
	}

bool IndexSet::
has_common_index(const IndexSet& other) const
    {
    for(int j = 1; j <= r_; ++j)
    for(int k = 1; k <= other.r_; ++k)
    if(index_[j] == other.index_[k]) return true;

    return false;
    }

bool IndexSet::
hasindex(const Index& I) const
	{
    if(I.m() == 1) return hasindex1(I);
    else           return hasindexn(I);
    return false;
	}

bool IndexSet::
hasindexn(const Index& I) const
	{
    for(int j = 1; j <= rn_; ++j)
    if(index_[j] == I) return true;
    return false;
	}

bool IndexSet::
hasindex1(const Index& I) const
	{
    for(int j = rn_+1; j <= r_; ++j)
    if(index_[j] == I) return true;
    return false;
	}

bool IndexSet::
hasAllIndex(const array<Index,NMAX+1>& I, int nind) const
    {
    for(int n = 1; n <= nind; ++n)
        {
        const Index& ii = I[n];
        bool found = false;
        if(ii.m() == 1)
            {
            for(int j = rn_+1; j <= r_; ++j)
                if(index_[j] == ii)
                    {
                    found = true;
                    break;
                    }
            }
        else
            {
            for(int j = 1; j <= rn_; ++j)
                if(index_[j] == ii)
                    {
                    found = true;
                    break;
                    }
            }
        if(!found) return false;
        }
    return true;
    }




void IndexSet::
getperm(const IndexSet& other, Permutation& P) const
	{
	for(int j = 1; j <= r_; ++j)
	    {
	    bool got_one = false;
	    for(int k = 1; k <= r_; ++k)
            if(other.index_[j] == index_[k])
                { P.from_to(j,k); got_one = true; break; }
	    if(!got_one)
            {
            std::cerr << "j = " << j << "\n";
            Print(*this); 
            std::cerr << "other.index_ = \n";
            for(int j = 1; j <= r_; ++j) 
                { std::cerr << other.index_[j] << "\n"; }
            Error("IndexSet::getperm: no matching index");
            }
	    }
	}

int IndexSet::
minM() const
    {
    if(rn_ == 0) return 1;

    int mm = index_[1].m();
    for(int j = 2; j <= rn_; ++j)
        mm = min(mm,index_[j].m());

    return mm;
    }

void IndexSet::
noprime(PrimeType p)
    {
    ur_ = 0;
    for(int j = 1; j <= r_; ++j) 
        {
        Index& J = index_[j];
        J.noprime(p);
        ur_ += J.uniqueReal();
        }
#ifdef SET_UR
        setUniqueReal();
#endif
	}

void IndexSet::
doprime(PrimeType pt, int inc)
	{
    ur_ = 0;
    for(int j = 1; j <= r_; ++j) 
        {
        Index& J = index_[j];
        J.doprime(pt,inc);
        ur_ += J.uniqueReal();
        }
#ifdef SET_UR
        setUniqueReal();
#endif
	}

void IndexSet::
mapprime(int plevold, int plevnew, PrimeType pt)
	{
    ur_ = 0;
    for(int j = 1; j <= r_; ++j) 
        {
        Index& J = index_[j];
        J.mapprime(plevold,plevnew,pt);
        ur_ += J.uniqueReal();
        }
#ifdef SET_UR
        setUniqueReal();
#endif
	}

void IndexSet::
mapprimeind(const Index& I, int plevold, int plevnew, PrimeType pt)
	{
    for(int j = (I.m() == 1 ? rn_+1 : 1); j <= r_; ++j) 
        if(index_[j] == I)
        {
        index_[j].mapprime(plevold,plevnew,pt);
        ur_ -= I.uniqueReal();
        ur_ += index_[j].uniqueReal();
#ifdef SET_UR
        setUniqueReal();
#endif
        return;
        }
    Print(*this);
    Print(I);
    Error("IndexSet::mapprimeind: index not found.");
	}

void IndexSet::
primeind(const Index& I, const Index& J)
	{ 
    mapindex(I,primed(I)); 
    mapindex(J,primed(J));
	}

/*
ITensor 
primeind(ITensor A, const Index& I1, const Index& I2)
    { 
    A.mapindex(I1,primed(I1));
    A.mapindex(I2,primed(I2));
    return A; 
    }
*/

//
// Methods for Manipulating IndexSets
//

void IndexSet::
mapindex(const Index& i1, const Index& i2)
	{
	assert(i1.m() == i2.m());
    if(i2.m() != i1.m())
        {
        Print(i1);
        Print(i2);
        Error("mapIndex: index must have matching m");
        }
	for(int j = 1; j <= r_; ++j) 
	    if(index_[j] == i1) 
		{
		index_[j] = i2;
        ur_ -= i1.uniqueReal();
        ur_ += i2.uniqueReal();
		return;
		}
	Print(i1);
	Error("IndexSet::mapindex: couldn't find i1.");
	}

void IndexSet::
addindex(const Index& I)
    {
    if(I.m() == 1)
        {
        index_[++r_] = I;
        }
    else
        {
#ifdef DEBUG
        if(r_ != rn_)
            Error("Adding m != 1 Index will overwrite m == 1 Index.");
#endif
        index_[++rn_] = I;
        ++r_;
        }
    ur_ += I.uniqueReal();
    }

void IndexSet::
addindexn(const array<Index,NMAX+1>& indices, int n) 
    {
#ifdef DEBUG
    if(r_ != rn_)
        Error("Adding m != 1 Index will overwrite m == 1 Index.");
#endif
    for(int j = 1; j <= n; ++j)
        {
        const Index& J = indices[j];
        index_[++rn_] = J;
        ur_ += J.uniqueReal();
        }
    r_ += n;
    }

void IndexSet::
addindexn(const Index& I)
    {
#ifdef DEBUG
    if(r_ != rn_)
        Error("Adding m != 1 Index will overwrite m == 1 Index.");
#endif
    index_[++rn_] = I;
    ur_ += I.uniqueReal();
    ++r_;
    }

void IndexSet::
addindex1(const array<Index,NMAX+1>& indices, int n) 
    {
    for(int j = 1; j <= n; ++j)
        {
        const Index& J = indices[j];
        index_[++r_] = J;
        ur_ += J.uniqueReal();
        }
    }

void IndexSet::
addindex1(const std::vector<Index>& indices) 
    { 
#ifdef DEBUG
    if((r_+(int)indices.size()) > NMAX)
        {
        Print(*this);
        Print(indices.size());
        Error("Too many indices added");
        }
#endif
    for(size_t j = 0; j < indices.size(); ++j)
        { 
        assert(indices[j].m() == 1);
        assert(!hasindex1(indices[j]));

        index_[++r_] = indices[j]; 
        ur_ += indices[j].uniqueReal();
        }
    }

void IndexSet::
addindex1(const Index& I) 
    { 
#ifdef DEBUG
    if(I.m() != 1)
        {
        Print(I);
        Error("Index must have m==1.");
        }
    if(hasindex1(I))
        {
        Print(*this);
        Print(I);
        Error("Adding Index twice");
        }
    if(r_ == NMAX) Error("Maximum number of indices reached");
#endif
    index_[++r_] = I;
    ur_ += I.uniqueReal();
    }

void IndexSet::
removeindex1(int j) 
    { 
    assert(j <= r_);
    assert(j > rn_);
    ur_ -= index_[j].uniqueReal();
    for(int k = j; k < r_; ++k) 
        index_[k] = index_[k+1];
    --r_;
    }

void IndexSet::
setUniqueReal()
	{
    ur_ = 0;
    for(int j = 1; j <= r_; ++j)
        ur_ += index_[j].uniqueReal();
	}

void IndexSet::
swap(IndexSet& other)
    {
    index_.swap(other.index_);

    int si = r_;
    r_ = other.r_;
    other.r_ = si;

    si = rn_;
    rn_ = other.rn_;
    other.rn_ = si;

    Real sr = ur_;
    ur_ = other.ur_;
    other.ur_ = sr;
    }

void IndexSet::
clear()
    {
    rn_ = 0;
    r_ = 0;
    ur_ = 0;
    }

std::ostream&
operator<<(std::ostream& s, const IndexSet& is)
    {
    int i = 1; 
    for(; i < is.r(); ++i) { s << is.index(i) << ", "; } 
    if(is.r() != 0) { s << is.index(i); } //print last one
    return s;
    }


void IndexSet::
read(std::istream& s)
    {
    s.read((char*) &r_,sizeof(r_));
    s.read((char*) &rn_,sizeof(rn_));
    ur_ = 0;
    for(int j = 1; j <= r_; ++j) 
        {
        index_[j].read(s);
        ur_ += index_[j].uniqueReal();
        }
    }

void IndexSet::
write(std::ostream& s) const
    {
    s.write((char*) &r_,sizeof(r_));
    s.write((char*) &rn_,sizeof(rn_));
    for(int j = 1; j <= r_; ++j) 
        index_[j].write(s);
    }
