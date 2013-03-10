//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEXSET_H
#define __ITENSOR_INDEXSET_H
#include "index.h"
#include "permutation.h"
#include "boost/make_shared.hpp"

#define Array boost::array
#define Cout std::cout
#define Endl std::endl

//
// IndexSet
//

template <class IndexT>
class IndexSet
    {
    public:

    IndexSet();

    explicit
    IndexSet(const IndexT& i1);

    IndexSet(const IndexT& i1, const IndexT& i2);

    IndexSet(IndexT i1, IndexT i2, IndexT i3,
             IndexT i4 = IndexT::Null(), 
             IndexT i5 = IndexT::Null(), 
             IndexT i6 = IndexT::Null(),
             IndexT i7 = IndexT::Null(), 
             IndexT i8 = IndexT::Null());

    template <class Iterable> 
    explicit
    IndexSet(const Iterable& ii, int size = -1, int offset = 0);

    template <class Iterable>
    IndexSet(const Iterable& ii, int size, int& alloc_size, int offset);

    IndexSet(const IndexSet& other, const Permutation& P);

    //
    // Type definitions
    //

    typedef Array<IndexT,NMAX>
    StorageType;

    typedef typename StorageType::const_iterator 
    const_iterator;

    typedef typename boost::shared_ptr<IndexSet<IndexT> >
    Ptr;

    //
    // Accessor Methods
    //

    int
    r() const { return r_; }

    int
    rn() const { return rn_; }

    const IndexT&
    index(int j) const { return index_[j-1]; }

    const IndexT&
    operator[](int j) const { return index_[j]; }

    int
    size() const { return r_; }

    int
    m(int j) const { return GET(index_,j-1).m(); }

    int
    dim() const;

    const_iterator
    begin() const { return index_.begin(); }

    const_iterator
    end() const { return (index_.begin()+r_); }

    const StorageType&
    storage() const { return index_; }

    Real
    uniqueReal() const { return ur_; }

    //
    // Index Analysis
    //

    bool 
    hastype(IndexType t) const;

    const IndexT&
    finddir(Arrow dir) const;

    int
    minM() const;

    int
    maxM() const;

    //
    // Primelevel Methods
    //

    void 
    prime(int inc = 1) { prime(All,inc); }

    void 
    prime(IndexType type, int inc = 1);

    void 
    prime(const IndexT& I, int inc = 1);

    void 
    noprime(IndexType type = All);

    void 
    noprime(const IndexT& I);

    void 
    mapprime(int plevold, int plevnew, IndexType type = All);

    void 
    mapprimeind(const IndexT& I, int plevold, int plevnew, 
                IndexType type = All);

    //
    // Methods for Manipulating IndexSets
    //

    void 
    addindex(const IndexT& I);

    void 
    addindexn(const Array<IndexT,NMAX+1>& indices, int n);

    void 
    addindexn(const IndexT& I);

    void 
    addindex1(const Array<IndexT,NMAX+1>& indices, int n);

    void 
    addindex1(const std::vector<IndexT>& indices);

    void 
    addindex1(const IndexT& I);

    void
    setUniqueReal();

    void
    swap(IndexSet& other);

    void
    clear();

    //
    // Other Methods
    //

    void
    reset();

    void
    conj();

    void
    read(std::istream& s);

    void
    write(std::ostream& s) const;

    static const Ptr& Null()
        {
        static Ptr Null_ = boost::make_shared<IndexSet<IndexT> >();
        return Null_;
        }

    private:

    //////////
    //
    // Data Members
    //

    StorageType index_;

    int rn_,
        r_;

    Real ur_;

    //
    /////////

    template <class Iterable>
    void
    sortIndices(const Iterable& I, int ninds, int& alloc_size, int offset = 0);

    };

template<class IndexT>
IndexSet<IndexT>::
IndexSet()
    :
    rn_(0),
    r_(0),
    ur_(0)
    { }

template<class IndexT>
IndexSet<IndexT>::
IndexSet(const IndexT& i1)
    :
    rn_((i1.m() == 1 ? 0 : 1)),
    r_(1),
    ur_(i1.uniqueReal())
    { 
#ifdef DEBUG
    if(i1 == IndexT::Null())
        Error("i1 is null");
#endif
    index_[0] = i1;
    }

template<class IndexT>
IndexSet<IndexT>::
IndexSet(const IndexT& i1, const IndexT& i2)
    :
    r_(2),
    ur_(i1.uniqueReal() + i2.uniqueReal())
    { 
#ifdef DEBUG
    if(i1 == IndexT::Null())
        Error("i1 is null");
    if(i2 == IndexT::Null())
        Error("i2 is null");
#endif
	if(i1.m()==1) 
	    {
	    index_[0] = i2; 
        index_[1] = i1; 
	    rn_ = (i2.m() == 1 ? 0 : 1);
	    }
	else 
	    { 
	    index_[0] = i1; 
        index_[1] = i2; 
	    rn_ = (i2.m() == 1 ? 1 : 2); 
	    }
    }

template<class IndexT>
IndexSet<IndexT>::
IndexSet(IndexT i1, IndexT i2, IndexT i3,
         IndexT i4, IndexT i5, IndexT i6,
         IndexT i7, IndexT i8)
    :
    r_(3)
    { 
#ifdef DEBUG
    if(i1 == IndexT::Null())
        Error("i1 is null");
    if(i2 == IndexT::Null())
        Error("i2 is null");
    if(i3 == IndexT::Null())
        Error("i3 is null");
#endif
	Array<IndexT,NMAX> ii = {{ i1, i2, i3, i4, i5, i6, i7, i8 }};
	while(ii[r_] != IndexT::Null()) ++r_;
    int alloc_size;
    sortIndices(ii,r_,alloc_size,0);
    setUniqueReal();
    }

template <class IndexT>
template <class Iterable>
IndexSet<IndexT>::
IndexSet(const Iterable& ii, int size, int offset)
    { 
    r_ = (size < 0 ? ii.size() : size);
    int alloc_size = -1;
    sortIndices(ii,r_,alloc_size,offset);
    setUniqueReal();
    }

template <class IndexT>
template <class Iterable>
IndexSet<IndexT>::
IndexSet(const Iterable& ii, int size, int& alloc_size, int offset)
    :
    r_(size)
    { 
    sortIndices(ii,size,alloc_size,offset);
    setUniqueReal();
    }


template <class IndexT>
IndexSet<IndexT>::
IndexSet(const IndexSet& other, const Permutation& P)
    :
    rn_(other.rn_),
    r_(other.r_),
    ur_(other.ur_)
    {
    for(int j = 1; j <= r_; ++j)
        index_[P.dest(j)-1] = other.index_[j-1];
    }

template <class IndexT>
int IndexSet<IndexT>::
dim() const
    {   
    int d = 1;
    for(int j = 0; j < rn_; ++j)
        d *= index_[j].m();
    return d;
    }

template <class IndexT>
bool IndexSet<IndexT>::
hastype(IndexType t) const
	{
    for(int j = 0; j < r_; ++j)
        if(index_[j].type() == t) return true;
    return false;
	}

template <class IndexT>
const IndexT& IndexSet<IndexT>::
finddir(Arrow dir) const
    {
    for(int j = 0; j < r_; ++j)
        if(index_[j].dir() == dir) return index_[j];
    Error("Couldn't find index with specified dir");
    return IndexT::Null();
    }



template <class IndexT>
int IndexSet<IndexT>::
minM() const
    {
    if(rn_ < r_) return 1;

    int mm = index_[0].m();
    for(int j = 1; j < rn_; ++j)
        mm = min(mm,index_[j].m());

    return mm;
    }

template <class IndexT>
int IndexSet<IndexT>::
maxM() const
    {
    if(rn_ == 0) return 1;

    int mm = index_[0].m();
    for(int j = 1; j < rn_; ++j)
        mm = max(mm,index_[j].m());

    return mm;
    }

template <class IndexT>
void IndexSet<IndexT>::
noprime(IndexType type)
    {
    ur_ = 0;
    for(int j = 0; j < r_; ++j) 
        {
        IndexT& J = index_[j];
#ifdef DEBUG
        //Check if calling noprime is ok
        //Error if it causes duplicate indices
        for(int k = 0; k < r_; ++k)
            {
            if(k != j && index_[j].noprimeEquals(index_[k]))
                {
                throw ITError("Calling noprime leads to duplicate indices");
                }
            }
#endif
        J.noprime(type);
        ur_ += J.uniqueReal();
        }
	}

template <class IndexT>
void IndexSet<IndexT>::
noprime(const IndexT& I)
    {
    int j = (I.m() == 1 ? rn_ : 0);
    for(; j < r_; ++j) 
        {
        if(index_[j] == I)
            {
#ifdef DEBUG
            //Check if calling noprime is ok
            //Error if it causes duplicate indices
            for(int k = 0; k < r_; ++k)
                {
                if(k != j && index_[j].noprimeEquals(index_[k]))
                    {
                    throw ITError("Calling noprime leads to duplicate indices");
                    }
                }
#endif
            index_[j].noprime();
            ur_ -= I.uniqueReal();
            ur_ += index_[j].uniqueReal();
            return;
            }
        }
    Print(*this);
    Print(I);
    Error("IndexSet::prime: index not found.");
    }

template <class IndexT>
void IndexSet<IndexT>::
prime(IndexType type, int inc)
	{
    ur_ = 0;
    for(int j = 0; j < r_; ++j) 
        {
        IndexT& J = index_[j];
        J.prime(type,inc);
        ur_ += J.uniqueReal();
        }
	}

template <class IndexT>
void IndexSet<IndexT>::
prime(const IndexT& I, int inc)
    {
    for(int j = (I.m() == 1 ? rn_ : 0); j < r_; ++j) 
        if(index_[j] == I)
        {
        index_[j].prime(inc);
        ur_ -= I.uniqueReal();
        ur_ += index_[j].uniqueReal();
        return;
        }
    Print(*this);
    Print(I);
    Error("IndexSet::prime: index not found.");
    }

template <class IndexT>
void IndexSet<IndexT>::
mapprime(int plevold, int plevnew, IndexType type)
	{
    ur_ = 0;
    for(int j = 0; j < r_; ++j) 
        {
        IndexT& J = index_[j];
        J.mapprime(plevold,plevnew,type);
        ur_ += J.uniqueReal();
        }
	}

template <class IndexT>
void IndexSet<IndexT>::
mapprimeind(const IndexT& I, int plevold, int plevnew, IndexType type)
	{
    for(int j = (I.m() == 1 ? rn_ : 0); j < r_; ++j) 
        if(index_[j] == I)
        {
        index_[j].mapprime(plevold,plevnew,type);
        ur_ -= I.uniqueReal();
        ur_ += index_[j].uniqueReal();
        return;
        }
    Print(*this);
    Print(I);
    Error("IndexSet::mapprimeind: index not found.");
	}


//
// Methods for Manipulating IndexSets
//

template <class IndexT>
void IndexSet<IndexT>::
addindex(const IndexT& I)
    {
#ifdef DEBUG
    if(I == IndexT::Null())
        Error("Index is null");
#endif
    if(I.m() == 1)
        addindex1(I);
    else
        addindexn(I);
    }

template <class IndexT>
void IndexSet<IndexT>::
addindexn(const Array<IndexT,NMAX+1>& indices, int n) 
    {
#ifdef DEBUG
    if(r_+n > NMAX) Error("Maximum number of indices reached");
#endif
    if(r_ != rn_)
        {
        //Move all m==1's over by n-1
        for(int k = r_+n-1; k >= rn_+n; --k)
            index_[k] = index_[k-1];
        }
    for(int j = 1; j <= n; ++j)
        {
        const IndexT& J = indices[j];
        index_[rn_] = J;
        ++rn_;
        ur_ += J.uniqueReal();
        }
    r_ += n;
    }

template <class IndexT>
void IndexSet<IndexT>::
addindexn(const IndexT& I)
    {
#ifdef DEBUG
    if(r_ == NMAX) Error("Maximum number of indices reached");
#endif
    if(r_ != rn_)
        {
        //Move all m==1's over by 1
        for(int k = r_; k > rn_; --k)
            index_[k] = index_[k-1];
        }
    index_[rn_] = I;
    ++rn_;
    ++r_;
    ur_ += I.uniqueReal();
    }

template <class IndexT>
void IndexSet<IndexT>::
addindex1(const Array<IndexT,NMAX+1>& indices, int n) 
    {
#ifdef DEBUG
    if(r_+n > NMAX) Error("Maximum number of indices reached");
#endif
    for(int j = 1; j <= n; ++j)
        {
        const IndexT& J = indices[j];
        index_[r_] = J;
        ++r_;
        ur_ += J.uniqueReal();
        }
    }

template <class IndexT>
void IndexSet<IndexT>::
addindex1(const std::vector<IndexT>& indices) 
    { 
#ifdef DEBUG
    if((r_+(int)indices.size()) > NMAX)
        {
        Print(*this);
        Print(indices.size());
        Error("Maximum number of indices reached");
        }
#endif
    for(size_t j = 0; j < indices.size(); ++j)
        { 
        assert(indices[j].m() == 1);

        index_[r_] = indices[j]; 
        ++r_;
        ur_ += indices[j].uniqueReal();
        }
    }

template <class IndexT>
void IndexSet<IndexT>::
addindex1(const IndexT& I) 
    { 
#ifdef DEBUG
    if(I.m() != 1)
        {
        Print(I);
        Error("Index must have m==1.");
        }
    for(int j = rn_; j < r_; ++j)
        if(index_[j] == I)
            {
            Print(*this);
            Print(I);
            Error("Adding Index twice");
            }
    if(r_ == NMAX) Error("Maximum number of indices reached");
#endif
    index_[r_] = I;
    ++r_;
    ur_ += I.uniqueReal();
    }

template <class IndexT>
void IndexSet<IndexT>::
setUniqueReal()
	{
    ur_ = 0;
    for(int j = 0; j < r_; ++j)
        ur_ += index_[j].uniqueReal();
	}

template <class IndexT>
void IndexSet<IndexT>::
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

template <class IndexT>
void IndexSet<IndexT>::
clear()
    {
    rn_ = 0;
    r_ = 0;
    ur_ = 0;
    }

template <class IndexT>
void IndexSet<IndexT>::
reset()
    {
    rn_ = 0;
    r_ = 0;
    ur_ = 0;
    }

template <class IndexT>
void IndexSet<IndexT>::
conj()
    {
    for(int j = 0; j < r_; ++j)
        index_[j].conj();
    }

template <class IndexT>
void IndexSet<IndexT>::
read(std::istream& s)
    {
    s.read((char*) &r_,sizeof(r_));
    s.read((char*) &rn_,sizeof(rn_));
    ur_ = 0;
    for(int j = 0; j < r_; ++j) 
        {
        index_[j].read(s);
        ur_ += index_[j].uniqueReal();
        }
    }

template <class IndexT>
void IndexSet<IndexT>::
write(std::ostream& s) const
    {
    s.write((char*) &r_,sizeof(r_));
    s.write((char*) &rn_,sizeof(rn_));
    for(int j = 0; j < r_; ++j) 
        index_[j].write(s);
    }

template <class IndexT>
template <class Iterable>
void IndexSet<IndexT>::
sortIndices(const Iterable& I, int ninds, int& alloc_size, int offset)
    {
#ifdef DEBUG
    if(ninds > NMAX)
        Error("Too many indices for IndexSet");
#endif

    rn_ = 0;
    alloc_size = 1;

    int r1_ = 0;
    Array<const IndexT*,NMAX> index1_;

    for(int n = offset; n < ninds+offset; ++n)
        {
        const IndexT& i = I[n];
#ifdef DEBUG
        if(i == IndexT::Null()) Error("Null Index in sortIndices");
#endif
        if(i.m()==1) 
            { 
            index1_[r1_] = &i;
            ++r1_;
            }
        else         
            { 
            index_[rn_] = i; 
            ++rn_;
            alloc_size *= i.m(); 
            }
        }
    for(int l = 0; l < r1_; ++l) 
        {
        index_[rn_+l] = *(index1_[l]);
        }
    }

template<class IndexT>
Arrow
dir(const IndexSet<IndexT>& is, const IndexT& I)
    {
    for(int j = 0; j < is.r(); ++j)
        {
        if(is[j] == I) 
            return is[j].dir();
        }
    Error("dir: Index not found");
    return In;
    }

template <class IndexT>
std::ostream&
operator<<(std::ostream& s, const IndexSet<IndexT>& is)
    {
    for(int i = 1; i <= is.r(); ++i) 
        s << is.index(i) << "\n"; 
    return s;
    }

template <> inline
std::ostream&
operator<<(std::ostream& s, const IndexSet<Index>& is)
    {
    int i = 1; 
    for(; i < is.r(); ++i) { s << is.index(i) << ", "; } 
    if(is.r() != 0) { s << is.index(i); } //print last one
    return s;
    }

//
// IndexSet helper methods
//

//
// Given IndexSet<IndexT> iset and IndexT I,
// return int j such that iset[j] == I.
// If not found, throws an ITError.
//
template <class IndexT>
int
findindex(const IndexSet<IndexT>& iset, const IndexT& I)
    {
    int j = (I.m()==1 ? iset.rn() : 0);
    for(; j < iset.r(); ++j)
        {
        if(iset[j] == I) return j;
        }
    Print(I);
    Error("Index I not found");
    return 0;
    }

template <class IndexT>
const IndexT&
findtype(const IndexSet<IndexT>& iset, IndexType t)
	{
    for(int j = 0; j < iset.r(); ++j)
        if(iset[j].type() == t) return iset[j];
    Error("findtype failed."); 
    return IndexT::Null();
	}

template <class IndexT, class ArrayType>
void
getperm(const IndexSet<IndexT>& iset, const ArrayType& ind, 
        Permutation& P)
	{
	for(int j = 0; j < iset.r(); ++j)
	    {
	    bool got_one = false;
	    for(int k = 0; k < iset.r(); ++k)
            {
            if(ind[j] == iset[k])
                { 
                P.from_to(j+1,k+1); 
                got_one = true; 
                break;
                }
            }
	    if(!got_one)
            {
            Cout << "j = " << j << "\n";
            Print(iset); 
            Cout << "ind = \n";
            for(int j = 0; j < iset.r(); ++j) 
                Cout << j << " " << ind[j] << "\n";
            Cout << Endl;
            Error("IndexSet::getperm: no matching index");
            }
	    }
	}

template <class IndexT>
bool
hasindex(const IndexSet<IndexT>& iset, const IndexT& I)
	{
    int j = (I.m()==1 ? iset.rn() : 0);
    for(; j < iset.r(); ++j)
        {
        if(iset[j] == I) return true;
        }
    return false;
	}


#undef Array
#undef Cout
#undef Endl

#endif
