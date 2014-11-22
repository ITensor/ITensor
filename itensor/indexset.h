//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEXSET_H
#define __ITENSOR_INDEXSET_H
#include "index.h"
#include "permutation.h"

namespace itensor {

template<typename IndexT>
class IndexSetT<IndexT>;

class IQIndex;

using IndexSet = IndexSetT<Index>;
using IQIndexSet = IndexSetT<IQIndex>;

//
// IndexSet
//

template <class IndexT>
class IndexSetT
    {
    public:

    using storage = std::vector<IndexT>;
    using const_iterator = typename storage::const_iterator;

    IndexSetT();

    explicit
    IndexSetT(const IndexT& i1);

    IndexSetT(const IndexT& i1, 
             const IndexT& i2);

    // construct from 3 or more Index's
    template <typename... Inds>
    IndexSetT(const IndexT& i1, 
             const Inds&... inds);

    IndexSetT(storage&& ii);

    template <class Iterable> 
    explicit
    IndexSetT(const Iterable& ii);

    template <class Iterable>
    IndexSetT(const Iterable& ii, size_t size, size_t offset);

    //
    // Accessor Methods
    //

    int
    r() const { return index_.size(); }
    
    size_t
    size() const { return index_.size(); }

    int
    rn() const { return rn_; }

    // 0-indexed access
    const IndexT&
    operator[](int j) const { return index_[j]; }

    // 1-indexed access
    const IndexT&
    index(int j) const { return index_[j-1]; }

    long
    dim() const;

    const IndexT&
    front() const { return index_.front(); }

    const IndexT&
    back() const { return index_.back(); }

    const_iterator
    begin() const { return index_.begin(); }

    const_iterator
    end() const { return index_.end(); }

    operator const storage&() const { return index_; }

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
    IndexSetT
    operator*(const IndexSetT& other) const;

    //
    // Other Methods
    //

    void 
    addindex(const IndexT& I);

    void
    replaceIndex(const IndexT& oind, const IndexT& nind);

    void
    swap(IndexSetT& other);

    void
    clear();

    void
    dag();

    void
    read(std::istream& s);

    void
    write(std::ostream& s) const;

    private:

    //////////

    storage index_;

    int rn_;

    /////////

    void
    init();

    bool static
    compare_index(const IndexT& i1,
                  const IndexT& i2)
        {
        return i1.m() > i2.m();
        }

    };

template<class IndexT>
IndexSetT<IndexT>::
IndexSetT()
    :
    rn_(0)
    { }

template<class IndexT>
IndexSetT<IndexT>::
IndexSetT(const IndexT& i1)
    :
    index_(1,i1),
    rn_((i1.m() == 1 ? 0 : 1))
    { 
#ifdef DEBUG
    if(i1 == IndexT::Null())
        Error("i1 is null");
#endif
    }

template<class IndexT>
IndexSetT<IndexT>::
IndexSetT(const IndexT& i1, 
         const IndexT& i2)
    :
    index_(2)
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
template<typename... Inds>
IndexSetT<IndexT>::
IndexSetT(const IndexT& i1, 
         const Inds&... inds)
    :
    index_{i1,inds...},
    rn_(0)
    { 
    init();
    }

template <class IndexT>
IndexSetT<IndexT>::
IndexSetT(storage&& ii)
    :
    index_(std::move(ii)),
    rn_(0)
    { 
    init();
    }

template <class IndexT>
template <class Iterable>
IndexSetT<IndexT>::
IndexSetT(const Iterable& ii)
    :
    index_(ii.begin(),ii.end()),
    rn_(0)
    { 
    init();
    }

template <class IndexT>
template <class Iterable>
IndexSetT<IndexT>::
IndexSetT(const Iterable& ii, size_t size, size_t offset)
    :
    index_(size),
    rn_(0)
    { 
    for(size_t n = offset, i = 0; n < size+offset; ++n, ++i)
        {
        index_[i] = ii[n];
        }
    init();
    }


template <class IndexT>
long IndexSetT<IndexT>::
dim() const
    {   
    long d = 1;
    for(int j = 0; j < rn_; ++j)
        d *= index_[j].m();
    return d;
    }

template <class IndexT>
void IndexSetT<IndexT>::
noprime(IndexType type)
    {
    for(auto& J : index_) J.noprime(type);

#ifdef DEBUG
        //Check if calling noprime is ok
        //Error if it causes duplicate indices
    for(size_t j = 0; j < index_.size(); ++j) 
        if(type == All || J.type() == type)
            {
            for(size_t k = 0; k < index_.size(); ++k)
                {
                const auto& K = index_[k];
                if(type != All && K.type() != type) continue;
                if(k != j && index_[j].noprimeEquals(index_[k]))
                    {
                    //Print(*this);
                    throw ITError("Calling noprime would lead to duplicate indices");
                    }
                }
            }
#endif
	}

template <class IndexT>
void IndexSetT<IndexT>::
noprime(const IndexT& I)
    {
    auto j = size_t(I.m() == 1 ? rn_ : 0);
    for(; j < index_.size(); ++j) 
        {
        if(index_[j] == I)
            {
#ifdef DEBUG
            //Check if calling noprime is ok
            //Error if it causes duplicate indices
            for(size_t k = 0; k < index_.size(); ++k)
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
    Error("IndexSetT::prime: index not found.");
    }

template <class IndexT>
void IndexSetT<IndexT>::
prime(IndexType type, int inc)
	{
    for(auto& J : index_) J.prime(type,inc);
	}

template <class IndexT>
void IndexSetT<IndexT>::
prime(const IndexT& I, 
      int inc)
    {
#ifdef DEBUG
    if(!I)
        {
        Error("Request to prime null index");
        }
#endif
    auto start = size_t(I.m() == 1 ? rn_ : 0);
    for(auto j = start; j < index_.size(); ++j) 
        if(index_[j] == I)
            {
            index_[j].prime(inc);
            return;
            }
    Print(*this);
    Print(I);
    Error("IndexSetT::prime: index not found.");
    }

template <class IndexT>
void IndexSetT<IndexT>::
mapprime(int plevold, int plevnew, IndexType type)
	{
    for(auto& J : index_) J.mapprime(plevold,plevnew,type);
	}

//template <class IndexT>
//IndexSetT<IndexT> inline IndexSetT<IndexT>::
//operator*(const IndexSetT& other) const
//    {
//    IndexSetT<IndexT> res;
//
//    //Loop over m!=1 indices of this
//    for(int i = 0; i < rn_; ++i)
//        {
//        const IndexT& I = index_[i];
//        //Loop over m!=1 indices of other
//        bool found = false;
//        for(int j = 0; j < other.rn_; ++j)
//            {
//            if(I == other.index_[j])
//                {
//                found = true;
//                break;
//                }
//            }
//        if(!found) 
//            res.addindex(I);
//        }
//
//    //Loop over m!=1 indices of other
//    for(int j = 0; j < other.rn_; ++j)
//        {
//        const IndexT& J = other.index_[j];
//        //Loop over m!=1 indices of other
//        bool found = false;
//        for(int i = 0; i < rn_; ++i)
//            {
//            if(J == index_[i])
//                {
//                found = true;
//                break;
//                }
//            }
//        if(!found) 
//            res.addindex(J);
//        }
//
//    //Loop over m==1 indices of this
//    for(int i = rn_; i < r_; ++i)
//        {
//        const IndexT& I = index_[i];
//        //Loop over m==1 indices of other
//        bool found = false;
//        for(int j = other.rn_; j < other.r_; ++j)
//            {
//            if(I == other.index_[j])
//                {
//                found = true;
//                break;
//                }
//            }
//        if(!found) 
//            res.addindex(I);
//        }
//
//    //Loop over m!=1 indices of other
//    for(int j = other.rn_; j < other.r_; ++j)
//        {
//        const IndexT& J = other.index_[j];
//        //Loop over m!=1 indices of other
//        bool found = false;
//        for(int i = rn_; i < r_; ++i)
//            {
//            if(J == index_[i])
//                {
//                found = true;
//                break;
//                }
//            }
//        if(!found) 
//            res.addindex(J);
//        }
//
//    return res;
//    }


//
// Methods for Manipulating IndexSetT
//

template <class IndexT>
void IndexSetT<IndexT>::
addindex(const IndexT& I)
    {
#ifdef DEBUG
    if(I == IndexT::Null()) Error("Index is null");

    for(const auto& J : index_)
        if(J == I)
            {
            Print(*this);
            Print(I);
            Error("Adding Index twice");
            }
#endif

    if(I.m() == 1) 
        {
        index_.push_back(I);
        }
    else
        {
        auto m_eq_one = [](const IndexT& ii) { return ii.m()==1; };
        auto it = std::find(index_.begin(),index_.end(),m_eq_one);
        index_.insert(it,I);
        ++rn_;
        }
    }

template <class IndexT>
void IndexSetT<IndexT>::
replaceIndex(const IndexT& oind, const IndexT& nind)
    {
    if(nind.m() != oind.m())
        {
        Print(nind);
        Print(oind);
        Error("replaceIndex: new index must have same dimension as old.");
        }
    bool found = false;
    for(auto& J : index_)
        {
        if(oind == J)
            {
            J = nind;
            found = true;
            break;
            }
        }
    if(!found)
        Error("replaceIndex: index not found");
    }

template <class IndexT>
void IndexSetT<IndexT>::
swap(IndexSetT& other)
    {
    index_.swap(other.index_);
    std::swap(rn_,other.rn_);
    }

template <class IndexT>
void IndexSetT<IndexT>::
clear()
    {
    rn_ = 0;
    index_.clear();
    }

template <class IndexT>
void IndexSetT<IndexT>::
dag()
    {
    for(auto& J : index_) J.dag();
    }

template <class IndexT>
void IndexSetT<IndexT>::
read(std::istream& s)
    {
    size_t size = 0;
    s.read((char*) &size,sizeof(size));
    index_.resize(size);
    for(auto& J : index_)
        {
        J.read(s);
        }

    s.read((char*) &rn_,sizeof(rn_));
    }

template <class IndexT>
void IndexSetT<IndexT>::
write(std::ostream& s) const
    {
    size_t size = index_.size();
    s.write((char*) &size,sizeof(size));
    for(auto& J : index_)
        {
        J.write(s);
        }

    s.write((char*) &rn_,sizeof(rn_));
    }


template<class IndexT>
IndexSetT<IndexT>::
init()
    {
#ifdef DEBUG
    if(rn_ != 0) Error("rn_ should be zero in init()");

    for(const auto& ii : index_)
        if(ii == IndexT::Null()) Error("Null index in IndexSetT");
#endif

    std::sort(index_.begin(),index_.end(),compare_index);

    for(const auto& ii : index_)
        {
        if(ii.m() == 1) break;
        ++rn_;
        }
    }

//template <class IndexT>
//template <class Iterable>
//void IndexSetT<IndexT>::
//sortIndices(const Iterable& I, int ninds, int& alloc_size, int offset)
//    {
//#ifdef DEBUG
//    if(ninds > NMAX)
//        Error("Too many indices for IndexSetT");
//#endif
//
//    rn_ = 0;
//    alloc_size = 1;
//
//    int r1_ = 0;
//    array<const IndexT*,NMAX> index1_;
//
//    for(int n = offset; n < ninds+offset; ++n)
//        {
//        const IndexT& i = I[n];
//#ifdef DEBUG
//        if(i == IndexT::Null()) Error("Null Index in sortIndices");
//#endif
//        if(i.m()==1) 
//            { 
//            index1_[r1_] = &i;
//            ++r1_;
//            }
//        else         
//            { 
//            index_[rn_] = i; 
//            ++rn_;
//            alloc_size *= i.m(); 
//            }
//        }
//    for(int l = 0; l < r1_; ++l) 
//        {
//        index_[rn_+l] = *(index1_[l]);
//        }
//    }

//
//
// IndexSetT helper methods
//
//

template<class IndexT>
Arrow
dir(const IndexSetT<IndexT>& is, const IndexT& I)
    {
    for(const auto& J : is)
        {
        if(J == I) return J.dir();
        }
    Error("dir: Index not found");
    return In;
    }


template <class IndexT>
const IndexT&
finddir(const IndexSetT<IndexT>& iset, Arrow dir)
    {
    for(const auto& J : iset)
        {
        if(J.dir() == dir) return J;
        }
    Error("Couldn't find index with specified dir");
    return IndexT::Null();
    }

//
// Given IndexSetT<IndexT> iset and IndexT I,
// return int j such that iset[j] == I.
// If not found, throws an ITError.
//
template <class IndexT>
int
findindex(const IndexSetT<IndexT>& iset, 
          const IndexT& I)
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
findtype(const IndexSetT<IndexT>& iset, IndexType t)
	{
    for(const auto& J : iset)
        {
        if(J.type() == t) return J;
        }
    Error("findtype failed."); 
    return IndexT::Null();
	}

//
// Compute the permutation P taking an IndexSetT iset
// to oset (of type IndexSetT or array<IndexT,NMAX>)
//
template <class IndexT>
void
getperm(const IndexSetT<IndexT>& iset, 
        const typename IndexSetT<IndexT>::storage& oset, 
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
            println("j = ",j);
            println("iset =");
            for(int j = 0; j < iset.r(); ++j)
                printfln("%d %s",j,iset[j]);
            println("\noset = ");
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
            throw ITError("IndexSetT::getperm: no matching index");
            }
	    }
	}

template <class IndexT>
bool
hasindex(const IndexSetT<IndexT>& iset, 
         const IndexT& I)
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
hastype(const IndexSetT<IndexT>& iset, 
        IndexType t)
	{
    for(const auto& J : iset)
        {
        if(J.type() == t) return true;
        }
    return false;
	}

template <class IndexT>
int
minM(const IndexSetT<IndexT>& iset)
    {
    if(iset.rn() < iset.r()) return 1;

    int mm = iset[0].m();
    for(int j = 1; j < iset.rn(); ++j)
        mm = min(mm,iset[j].m());

    return mm;
    }

template <class IndexT>
int
maxM(const IndexSetT<IndexT>& iset)
    {
    if(iset.rn() == 0) return 1;

    int mm = iset[0].m();
    for(int j = 1; j < iset.rn(); ++j)
        mm = max(mm,iset[j].m());

    return mm;
    }

template <class IndexT>
std::ostream&
operator<<(std::ostream& s, const IndexSetT<IndexT>& is)
    {
    for(const auto& J : iset)
        {
        s << J << "\n"; 
        }
    return s;
    }

template <> inline
std::ostream&
operator<<(std::ostream& s, const IndexSetT<Index>& is)
    {
    int i = 1; 
    for(; i < is.r(); ++i) { s << is.index(i) << ", "; } 
    if(is.r() != 0) { s << is.index(i); } //print last one
    return s;
    }

}; //namespace itensor

#endif
