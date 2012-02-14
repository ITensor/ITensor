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
    r_(1)
    { 
    index_[1] = i1;
    setUniqueReal();
    }

IndexSet::
IndexSet(const Index& i1, const Index& i2)
    :
    r_(2)
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
    setUniqueReal();
    }

void IndexSet::
setUniqueReal()
	{
    ur_ = 0;
    for(int j = 1; j <= r_; ++j)
        ur_ += index_[j].uniqueReal();
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
        }
    setUniqueReal();
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
    setUniqueReal();
    }

void IndexSet::
removeindex1(int j) 
    { 
    assert(j <= r_);
    assert(j > rn_);
    for(int k = j; k < r_; ++k) 
        index_[k] = index_[k+1];
    --r_;
    setUniqueReal();
    }

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
	    if(GET(index_,j) == i1) 
		{
		GET(index_,j) = i2;
		setUniqueReal();
		return;
		}
	Print(i1);
	Error("IndexSet::mapindex: couldn't find i1.");
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

void IndexSet::
noprime(PrimeType p)
    {
    for(int j = 1; j <= r_; ++j) 
        index_[j].noprime(p);
    setUniqueReal();
	}

void IndexSet::
doprime(PrimeType pt, int inc)
	{
    for(int j = 1; j <= r_; ++j) 
        index_[j].doprime(pt,inc);
    setUniqueReal();
	}

void IndexSet::
mapprime(int plevold, int plevnew, PrimeType pt)
	{
    for(int j = 1; j <= r_; ++j) 
        index_[j].mapprime(plevold,plevnew,pt);
    setUniqueReal();
	}

void IndexSet::
mapprimeind(const Index& I, int plevold, int plevnew, PrimeType pt)
	{
    for(int j = (I.m() == 1 ? rn_+1 : 1); j <= r_; ++j) 
        if(index_[j] == I)
        {
        index_[j].mapprime(plevold,plevnew,pt);
        setUniqueReal();
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

std::ostream&
operator<<(std::ostream& s, const IndexSet& is)
    {
    int i = 1; 
    for(; i < is.r(); ++i) { s << is.index(i) << ", "; } 
    if(is.r() != 0) { s << is.index(i); } //print last one
    return s;
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
read(std::istream& s)
    {
    s.read((char*) &r_,sizeof(r_));
    s.read((char*) &rn_,sizeof(rn_));
    for(int j = 1; j <= r_; ++j) 
        index_[j].read(s);
    setUniqueReal();
    }

void IndexSet::
write(std::ostream& s) const
    {
    s.write((char*) &r_,sizeof(r_));
    s.write((char*) &rn_,sizeof(rn_));
    for(int j = 1; j <= r_; ++j) 
        index_[j].write(s);
    }
