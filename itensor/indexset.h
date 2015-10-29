//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEXSET_H
#define __ITENSOR_INDEXSET_H
#include "itensor/index.h"
#include "itensor/permutation.h"

namespace itensor {


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

    using Storage = array<IndexT,NMAX>;

    using const_iterator = typename Storage::const_iterator;

    using Ptr = shared_ptr<IndexSet<IndexT>>;

    //
    // Accessor Methods
    //

    int
    r() const { return r_; }

    int
    rn() const { return rn_; }

    int
    size() const { return r_; }

    const IndexT&
    index(int j) const { return index_[j-1]; }

    const IndexT&
    operator[](int j) const { return index_[j]; }

    int
    dim() const;

    IndexT
    front() const;

    IndexT
    back() const;

    const_iterator
    begin() const { return index_.begin(); }

    const_iterator
    end() const { return (index_.begin()+r_); }

    operator const Storage&() const { return index_; }

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

    //
    // Operators
    //

    //Contraction - just like tensor contraction but only the indices,
    //no data involved. Result is disjoint union of this and other
    //(this U other - this N other, where N is intersection).
    IndexSet
    operator*(const IndexSet& other) const;

    //
    // Other Methods
    //

    void 
    addindex(const IndexT& I);

    void
    replaceIndex(const IndexT& oind, const IndexT& nind);

    void
    swap(IndexSet& other);

    void
    clear();

    void
    dag();

    void
    read(std::istream& s);

    void
    write(std::ostream& s) const;

    //static const Ptr& Null()
    //    {
    //    static Ptr Null_ = make_shared<IndexSet<IndexT> >();
    //    return Null_;
    //    }

    private:

    //////////

    Storage index_;

    int rn_,
        r_;

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
    r_(0)
    { }

template<class IndexT>
IndexSet<IndexT>::
IndexSet(const IndexT& i1)
    :
    rn_((i1.m() == 1 ? 0 : 1)),
    r_(1)
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
    r_(2)
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
    array<IndexT,NMAX> ii = {{ i1, i2, i3, i4, i5, i6, i7, i8 }};
	while(r_ < NMAX && ii[r_] != IndexT::Null()) ++r_;
    int alloc_size;
    sortIndices(ii,r_,alloc_size,0);
    }

template <class IndexT>
template <class Iterable>
IndexSet<IndexT>::
IndexSet(const Iterable& ii, int size, int offset)
    { 
    r_ = (size < 0 ? ii.size() : size);
    int alloc_size = -1;
    sortIndices(ii,r_,alloc_size,offset);
    }

template <class IndexT>
template <class Iterable>
IndexSet<IndexT>::
IndexSet(const Iterable& ii, int size, int& alloc_size, int offset)
    :
    r_(size)
    { 
    sortIndices(ii,size,alloc_size,offset);
    }


template <class IndexT>
IndexSet<IndexT>::
IndexSet(const IndexSet& other, const Permutation& P)
    :
    rn_(other.rn_),
    r_(other.r_)
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
IndexT IndexSet<IndexT>::
front() const
    {
#ifdef DEBUG
    if(r_ == 0)
        Error("Empty IndexSet");
#endif
    return index_[0];
    }

template <class IndexT>
IndexT IndexSet<IndexT>::
back() const
    {
#ifdef DEBUG
    if(r_ == 0)
        Error("Empty IndexSet");
#endif
    return index_[r_-1];
    }

template <class IndexT>
void IndexSet<IndexT>::
noprime(IndexType type)
    {
    for(int j = 0; j < r_; ++j) 
        {
        IndexT& J = index_[j];
#ifdef DEBUG
        //Check if calling noprime is ok
        //Error if it causes duplicate indices
        if(type == All || J.type() == type)
            {
            for(int k = 0; k < r_; ++k)
                {
                const IndexT& K = index_[k];
                if(type != All && K.type() != type) continue;
                if(k != j && index_[j].noprimeEquals(index_[k]))
                    {
                    //Print(*this);
                    throw ITError("Calling noprime would lead to duplicate indices");
                    }
                }
            }
#endif
        J.noprime(type);
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
    for(int j = 0; j < r_; ++j) 
        {
        IndexT& J = index_[j];
        J.prime(type,inc);
        }
	}

template <class IndexT>
void IndexSet<IndexT>::
prime(const IndexT& I, int inc)
    {
#ifdef DEBUG
    if(!I)
        {
        Error("Request to prime null index");
        }
#endif
    for(int j = (I.m() == 1 ? rn_ : 0); j < r_; ++j) 
        if(index_[j] == I)
        {
        index_[j].prime(inc);
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
    for(int j = 0; j < r_; ++j) 
        {
        IndexT& J = index_[j];
        J.mapprime(plevold,plevnew,type);
        }
	}

template <class IndexT>
IndexSet<IndexT> inline IndexSet<IndexT>::
operator*(const IndexSet& other) const
    {
    IndexSet<IndexT> res;

    //Loop over m!=1 indices of this
    for(int i = 0; i < rn_; ++i)
        {
        const IndexT& I = index_[i];
        //Loop over m!=1 indices of other
        bool found = false;
        for(int j = 0; j < other.rn_; ++j)
            {
            if(I == other.index_[j])
                {
                found = true;
                break;
                }
            }
        if(!found) 
            res.addindex(I);
        }

    //Loop over m!=1 indices of other
    for(int j = 0; j < other.rn_; ++j)
        {
        const IndexT& J = other.index_[j];
        //Loop over m!=1 indices of other
        bool found = false;
        for(int i = 0; i < rn_; ++i)
            {
            if(J == index_[i])
                {
                found = true;
                break;
                }
            }
        if(!found) 
            res.addindex(J);
        }

    //Loop over m==1 indices of this
    for(int i = rn_; i < r_; ++i)
        {
        const IndexT& I = index_[i];
        //Loop over m==1 indices of other
        bool found = false;
        for(int j = other.rn_; j < other.r_; ++j)
            {
            if(I == other.index_[j])
                {
                found = true;
                break;
                }
            }
        if(!found) 
            res.addindex(I);
        }

    //Loop over m!=1 indices of other
    for(int j = other.rn_; j < other.r_; ++j)
        {
        const IndexT& J = other.index_[j];
        //Loop over m!=1 indices of other
        bool found = false;
        for(int i = rn_; i < r_; ++i)
            {
            if(J == index_[i])
                {
                found = true;
                break;
                }
            }
        if(!found) 
            res.addindex(J);
        }

    return res;
    }


//
// Methods for Manipulating IndexSets
//

template <class IndexT>
void IndexSet<IndexT>::
addindex(const IndexT& I)
    {
#ifdef DEBUG
    if(r_ == NMAX) 
        Error("Maximum number of indices reached");
    if(I == IndexT::Null())
        Error("Index is null");
    for(int j = (I.m()==1 ? rn_ : 0); j < r_; ++j)
        {
        if(index_[j] == I)
            {
            Print(*this);
            Print(I);
            Error("Adding Index twice");
            }
        }
#endif
    if(I.m() == 1)
        {
        index_[r_] = I;
        }
    else
        {
        if(r_ != rn_)
            {
            //Move all m==1's over by 1
            for(int k = r_; k > rn_; --k)
                index_[k] = index_[k-1];
            }
        index_[rn_] = I;
        ++rn_;
        }
    ++r_;
    }

template <class IndexT>
void IndexSet<IndexT>::
replaceIndex(const IndexT& oind, const IndexT& nind)
    {
    if(nind.m() != oind.m())
        {
        Print(nind);
        Print(oind);
        Error("replaceIndex: new index must have same dimension as old.");
        }
    bool found = false;
    for(int j = 0; j < r_; ++j) 
        {
        if(index_[j] == oind)
            {
            index_[j] = nind;
            found = true;
            }
        }
    if(!found)
        Error("replaceIndex: index not found");
    }

/*
template <class IndexT>
void IndexSet<IndexT>::
addindex1(const array<IndexT,NMAX+1>& indices, int n) 
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
    */

/*
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
    */

template <class IndexT>
void IndexSet<IndexT>::
swap(IndexSet& other)
    {
    index_.swap(other.index_);

    int tmp = r_;
    r_ = other.r_;
    other.r_ = tmp;

    tmp = rn_;
    rn_ = other.rn_;
    other.rn_ = tmp;
    }

template <class IndexT>
void IndexSet<IndexT>::
clear()
    {
    rn_ = 0;
    r_ = 0;
    }

template <class IndexT>
void IndexSet<IndexT>::
dag()
    {
    for(int j = 0; j < r_; ++j)
        index_[j].dag();
    }

template <class IndexT>
void IndexSet<IndexT>::
read(std::istream& s)
    {
    s.read((char*) &r_,sizeof(r_));
    s.read((char*) &rn_,sizeof(rn_));
    for(int j = 0; j < r_; ++j) 
        {
        index_[j].read(s);
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
    array<const IndexT*,NMAX> index1_;

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

//
//
// IndexSet helper methods
//
//

template<class IndexT>
int
rank(const IndexSet<IndexT>& is) { return is.r(); }

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
const IndexT&
finddir(const IndexSet<IndexT>& iset, Arrow dir)
    {
    for(int j = 0; j < iset.r(); ++j)
        if(iset[j].dir() == dir) return iset[j];
    Error("Couldn't find index with specified dir");
    return IndexT::Null();
    }

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

//
// Compute the permutation P taking an IndexSet iset
// to oset (of type IndexSet or array<IndexT,NMAX>)
//
template <class IndexT>
void
getperm(const IndexSet<IndexT>& iset, 
        const typename IndexSet<IndexT>::Storage& oset, 
        Permutation& P)
	{
	for(int j = 0; j < iset.r(); ++j)
	    {
	    bool got_one = false;
	    for(int k = 0; k < iset.r(); ++k)
            {
            if(oset[j] == iset[k])
                { 
                P.setFromTo(j+1,k+1); 
                got_one = true; 
                break;
                }
            }
	    if(!got_one)
            {
            println("Index sets are not permutations of each other");
            println("j = ",j);
            println("index set 1 =");
            for(int j = 0; j < iset.r(); ++j)
                printfln("%d %s",j,iset[j]);
            println("\nindex set 2 = ");
            for(int j = 0; j < iset.r(); ++j)
                printfln("%d %s",j,oset[j]);
            println();
            //printfln("iset uniqueReal = %.15E",iset.uniqueReal());
            //Real our = 0;
            //for(int i = 0; i < iset.r(); ++i)
            //    {
            //    our += oset[i].uniqueReal();
            //    }
            //printfln("oset uniqueReal = %.15E",our);
            //printfln("uniqueReal diff = %.15E",fabs(our-iset.uniqueReal()));
            throw ITError("IndexSet::getperm: no matching index");
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

template <class IndexT>
bool
hastype(const IndexSet<IndexT>& iset, IndexType t)
	{
    for(int j = 0; j < iset.r(); ++j)
        if(iset[j].type() == t) return true;
    return false;
	}

template <class IndexT>
int
minM(const IndexSet<IndexT>& iset)
    {
    if(iset.rn() < iset.r()) return 1;

    int mm = iset[0].m();
    for(int j = 1; j < iset.rn(); ++j)
        mm = min(mm,iset[j].m());

    return mm;
    }

template <class IndexT>
int
maxM(const IndexSet<IndexT>& iset)
    {
    if(iset.rn() == 0) return 1;

    int mm = iset[0].m();
    for(int j = 1; j < iset.rn(); ++j)
        mm = max(mm,iset[j].m());

    return mm;
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
    for(; i < is.r(); ++i) { s << is.index(i) << " "; } 
    if(is.r() != 0) { s << is.index(i); } //print last one
    return s;
    }

} //namespace itensor

#endif
