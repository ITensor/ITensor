//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "iqindexset.h"
using namespace std;
using boost::format;
using boost::array;

void 
intrusive_ptr_add_ref(IQIndexSet* p) 
    { 
    ++(p->numref); 
    }

void 
intrusive_ptr_release(IQIndexSet* p) 
    { 
    if(--(p->numref) == 0)
        { 
        delete p; 
        } 
    }

IQIndexSet::
IQIndexSet()
    :
    ur_(0),
    numref(0)
    { }

IQIndexSet::
IQIndexSet(const IQIndex& i1)
    :
    index_(1),
    ur_(i1.uniqueReal()),
    numref(0)
    { 
#ifdef DEBUG
    if(i1 == IQIndex::Null())
        Error("i1 was null");
#endif
    index_[0] = i1;
    }

IQIndexSet::
IQIndexSet(const IQIndex& i1, const IQIndex& i2)
    :
    index_(2),
    ur_(i1.uniqueReal() + i2.uniqueReal()),
    numref(0)
    { 
#ifdef DEBUG
    if(i1 == IQIndex::Null())
        Error("i1 was null");
    if(i2 == IQIndex::Null())
        Error("i2 was null");
#endif
    index_[0] = i1;
    index_[1] = i2;
    }

IQIndexSet::
IQIndexSet(const IQIndex& i1, const IQIndex& i2,
           const IQIndex& i3)
    :
    index_(3),
    ur_(i1.uniqueReal() + i2.uniqueReal() + i3.uniqueReal()),
    numref(0)
    { 
#ifdef DEBUG
    if(i1 == IQIndex::Null())
        Error("i1 was null");
    if(i2 == IQIndex::Null())
        Error("i2 was null");
    if(i3 == IQIndex::Null())
        Error("i3 was null");
#endif
    index_[0] = i1;
    index_[1] = i2;
    index_[2] = i3;
    }

IQIndexSet::
IQIndexSet(IQIndex i1, IQIndex i2, IQIndex i3,
         IQIndex i4, IQIndex i5, IQIndex i6,
         IQIndex i7, IQIndex i8)
    :
    index_(4),
    numref(0)
    { 
#ifdef DEBUG
    if(i1 == IQIndex::Null())
        Error("i1 was null");
    if(i2 == IQIndex::Null())
        Error("i2 was null");
    if(i3 == IQIndex::Null())
        Error("i3 was null");
    if(i4 == IQIndex::Null())
        Error("i4 was null");
#endif
	index_[0] = i1; 
	index_[1] = i2; 
	index_[2] = i3; 
	index_[3] = i4; 

	if(i5 != IQIndex::Null()) 
	    index_.push_back(i5);
	if(i6 != IQIndex::Null()) 
	    index_.push_back(i6);
	if(i7 != IQIndex::Null()) 
	    index_.push_back(i7);
	if(i8 != IQIndex::Null()) 
	    index_.push_back(i8);

    setUniqueReal();
    }

IQIndexSet::
IQIndexSet(std::vector<IQIndex>& iqinds)
    :
    numref(0)
    {
    index_.swap(iqinds);
    setUniqueReal();
    }

IQIndexSet::
IQIndexSet(const IQIndexSet& other)
    :
    index_(other.index_),
    ur_(other.ur_),
    numref(0)
    {
    }

IQIndexSet::
IQIndexSet(std::istream& s) 
    : 
    numref(0)
    { 
    int size;
    s.read((char*) &size,sizeof(size));
    index_.resize(size);
    ur_ = 0;
    for(int j = 0; j < size; ++j) 
        {
        index_[j].read(s);
        ur_ += index_[j].uniqueReal();
        }
    }

IQIndexSet::
IQIndexSet(int init_numref)
    :
    ur_(0),
    numref(init_numref)
    { }

const IQIndex& IQIndexSet::
findtype(IndexType t) const
	{
    for(size_t j = 0; j < index_.size(); ++j)
        if(index_[j].type() == t) return index_[j];
    Error("IQIndexSet::findtype failed."); 
    return IQIndex::Null();
	}

bool IQIndexSet::
findtype(IndexType t, IQIndex& I) const
	{
    for(size_t j = 0; j < index_.size(); ++j)
    if(index_[j].type() == t)
        {
        I = index_[j];
        return true;
        }
    return false;
	}

int IQIndexSet::
findindex(const IQIndex& I) const
    {
    for(int j = 0; j < (int)index_.size(); ++j)
        if(index_[j] == I) return (j+1);
    return 0;
	}

const IQIndex& IQIndexSet::
finddir(Arrow dir) const
    {
    Foreach(const IQIndex& I, index_)
        if(I.dir() == dir) return I;
    Error("Couldn't find arrow with specified dir");
    return IQIndex::Null();
    }

bool IQIndexSet::
has_common_index(const IQIndexSet& other) const
    {
    for(size_t j = 0; j < index_.size(); ++j)
    for(size_t k = 0; k <= other.index_.size(); ++k)
    if(index_[j] == other.index_[k]) return true;

    return false;
    }

bool IQIndexSet::
hasindex(const IQIndex& I) const
	{
    for(size_t j = 0; j < index_.size(); ++j)
    if(index_[j] == I) return true;
    return false;
	}

bool IQIndexSet::
hastype(IndexType t) const
    {
    for(size_t j = 0; j < index_.size(); ++j)
    if(index_[j].type() == t) return true;
    return false;
    }

int IQIndexSet::
minM() const
    {
    if(index_.empty()) return 1;

    int mm = index_[0].m();
    for(size_t j = 1; j < index_.size(); ++j)
        mm = min(mm,index_[j].m());

    return mm;
    }

int IQIndexSet::
maxM() const
    {
    if(index_.empty()) return 1;

    int mm = index_[0].m();
    for(size_t j = 1; j < index_.size(); ++j)
        mm = max(mm,index_[j].m());

    return mm;
    }

void IQIndexSet::
noprime(PrimeType p)
    {
    ur_ = 0;
    for(size_t j = 0; j < index_.size(); ++j) 
        {
        IQIndex& J = index_[j];
        J.noprime(p);
        ur_ += J.uniqueReal();
        }
	}

void IQIndexSet::
doprime(PrimeType pt, int inc)
	{
    ur_ = 0;
    for(size_t j = 0; j < index_.size(); ++j) 
        {
        IQIndex& J = index_[j];
        J.doprime(pt,inc);
        ur_ += J.uniqueReal();
        }
	}

void IQIndexSet::
mapprime(int plevold, int plevnew, PrimeType pt)
	{
    ur_ = 0;
    for(size_t j = 0; j < index_.size(); ++j) 
        {
        IQIndex& J = index_[j];
        J.mapprime(plevold,plevnew,pt);
        ur_ += J.uniqueReal();
        }
	}

void IQIndexSet::
mapprimeind(const IQIndex& I, int plevold, int plevnew, PrimeType pt)
	{
    for(size_t j = 0; j < index_.size(); ++j) 
        if(index_[j] == I)
        {
        index_[j].mapprime(plevold,plevnew,pt);
        ur_ -= I.uniqueReal();
        ur_ += index_[j].uniqueReal();
        return;
        }
    //Print(*this);
    //Print(I);
    Error("IQIndexSet::mapprimeind: index not found.");
	}

void IQIndexSet::
primeind(const IQIndex& I, const IQIndex& J)
	{ 
    mapindex(I,primed(I)); 
    mapindex(J,primed(J));
	}

void IQIndexSet::
indIncPrime(const IQIndex& I, int inc)
    {
    Foreach(IQIndex& J, index_)
        if(J.noprime_equals(I))
        {
        int p = J.primeLevel();
        J.mapprime(p,p+inc);
        return;
        }
    cerr << "IQIndex was " << I << "\n";
    cout << "IQIndex was " << I << "\n";
    Error("indIncPrime: couldn't find IQIndex");
    }

void IQIndexSet::
noprimeind(const IQIndex& I)
    {
#ifdef DEBUG
    int nmatch = 0;
	for(size_t j = 0; j < index_.size(); ++j) 
	    if(index_[j].noprime_equals(I))
            {
            ++nmatch;
            }
    if(nmatch > 1)
        {
        Print(*this);
        Print(I);
        Error("Calling noprimeind would lead to two copies of I");
        }
#endif

	for(size_t j = 0; j < index_.size(); ++j) 
	    if(index_[j] == I) 
		{
        ur_ -= index_[j].uniqueReal();
		index_[j].noprime();
        ur_ += index_[j].uniqueReal();
		return;
		}
    Print(I);
    Error("IQIndexSet does not contain IQIndex I");
    }

//
// Methods for Manipulating IQIndexSets
//

void IQIndexSet::
mapindex(const IQIndex& i1, const IQIndex& i2)
	{
    if(i2.m() != i1.m())
        {
        Print(i1);
        Print(i2);
        Error("mapIndex: index must have matching m");
        }
	for(size_t j = 0; j < index_.size(); ++j) 
	    if(index_[j] == i1) 
		{
		index_[j] = i2;
        ur_ -= i1.uniqueReal();
        ur_ += i2.uniqueReal();
		return;
		}
	Print(i1);
	Error("IQIndexSet::mapindex: couldn't find i1.");
	}

void IQIndexSet::
addindex(const IQIndex& I)
    {
#ifdef DEBUG
    if(I == IQIndex::Null())
        Error("IQIndex is null");
#endif
    index_.push_back(I);
    ur_ += I.uniqueReal();
    }

void IQIndexSet::
removeindex(int j) 
    { 
    ur_ -= index_[j].uniqueReal();
    index_.erase(index_.begin() + (j-1));
    }

void IQIndexSet::
setUniqueReal()
	{
    ur_ = 0;
    for(size_t j = 0; j < index_.size(); ++j)
        ur_ += index_[j].uniqueReal();
	}

void IQIndexSet::
swap(IQIndexSet& other)
    {
    index_.swap(other.index_);

    Real sr = ur_;
    ur_ = other.ur_;
    other.ur_ = sr;
    }

void IQIndexSet::
swapInds(vector<IQIndex>& newinds)
    {
    index_.swap(newinds);
    setUniqueReal();
    }

void IQIndexSet::
clear()
    {
    index_.clear();
    ur_ = 0;
    }

std::ostream&
operator<<(std::ostream& s, const IQIndexSet& is)
    {
    for(int i = 1; i <= is.r(); ++i) 
        s << is.index(i) << "\n"; 
    return s;
    }

void IQIndexSet::
conj()
    {
    for(size_t j = 0; j < index_.size(); ++j)
        index_[j].conj();
    }

void IQIndexSet::
conj(const IQIndex& I)
    {
    for(size_t j = 0; j < index_.size(); ++j)
        if(index_[j] == I)
            {
            index_[j].conj();
            return;
            }
    }

void IQIndexSet::
write(std::ostream& s) const
    {
    int size = index_.size();
    s.write((char*) &size,sizeof(size));
    for(int j = 0; j < size; ++j) 
        index_[j].write(s);
    }
