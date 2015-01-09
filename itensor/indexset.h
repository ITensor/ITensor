//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEXSET_H
#define __ITENSOR_INDEXSET_H
#include "index.h"
#include "permutation.h"
#include "range.h"

namespace itensor {

template <class IndexT>
class IndexSetT;

class IQIndex;

//
// IndexSetT
//
// Aliases:
using IndexSet = IndexSetT<Index>;
using IQIndexSet = IndexSetT<IQIndex>;

template <class IndexT>
class IndexSetT
    {
    public:

    using storage = std::vector<IndexT>;
    using const_iterator = typename storage::const_iterator;

    IndexSetT();

    explicit
    IndexSetT(const IndexT& i1);

    // construct from 2 or more indices
    template <typename... Inds>
    IndexSetT(const IndexT& i1, 
              const IndexT& i2,
              const Inds&... inds);

    template <class Iterable> 
    explicit
    IndexSetT(const Iterable& ii);

    template <class Iterable>
    IndexSetT(const Iterable& ii, size_t size, size_t offset);

    //
    // Accessor Methods
    //

    long
    dim(long i) const { return index_[i].m(); }
    long
    stride(long i) const { return stride_[i]; }

    int
    r() const { return index_.size(); }
    
    size_t
    size() const { return index_.size(); }

    bool
    empty() const { return index_.empty(); }

    int
    rn() const { return rn_; }

    // 0-indexed access
    const IndexT&
    operator[](int j) const { return index_[j]; }

    // 1-indexed access
    const IndexT&
    index(int j) const { return index_[j-1]; }

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
    // Other Methods
    //

    bool
    operator==(const IndexSetT& other) const;
    bool
    operator!=(const IndexSetT& other) const { return !operator==(other); }

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
    std::vector<long> stride_;
    int rn_;

    /////////

    void
    init();

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
    stride_(1,1),
    rn_((i1.m() == 1 ? 0 : 1))
    { 
#ifdef DEBUG
    if(!i1) Error("i1 is default initialized");
#endif
    }

template<class IndexT>
template<typename... Inds>
IndexSetT<IndexT>::
IndexSetT(const IndexT& i1, 
         const IndexT& i2,
         const Inds&... inds)
    :
    index_{i1,i2,inds...},
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
void IndexSetT<IndexT>::
noprime(IndexType type)
    {
    for(auto& J : index_) J.noprime(type);

#ifdef DEBUG
        //Check if calling noprime is ok
        //Error if it causes duplicate indices
    for(size_t j = 0; j < index_.size(); ++j) 
        if(type == All || index_[j].type() == type)
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

template <class IndexT>
bool IndexSetT<IndexT>::
operator==(const IndexSetT& other) const
    {
    if(other.r() != r()) return false;

    //IndexSetT sorts its indices by dimension
    //and id number, so enough to check if exactly
    //the same to check if indices are the same
    //in an unordered sense
    for(size_t j = 0; j < index_.size(); ++j)
        {
        if(index_[j] != other.index_[j]) return false;
        }
    return true;
    }


//
// Methods for Manipulating IndexSetT
//

template <class IndexT>
void IndexSetT<IndexT>::
addindex(const IndexT& I)
    {
#ifdef DEBUG
    if(!I) Error("Index is default initialized");

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
void IndexSetT<IndexT>::
init()
    {
#ifdef DEBUG
    if(rn_ != 0) Error("rn_ should be zero in init()");

    for(const auto& ii : index_)
        if(!ii) Error("Default initialized index in IndexSetT");
#endif

    auto comp = [](const IndexT& i1, const IndexT& i2) { return i1 > i2; };
    std::sort(index_.begin(),index_.end(),comp);

    stride_ = std::vector<long>(index_.size(),1);

    long str = 1;
    for(size_t j = 0; j < index_.size(); ++j)
        {
        stride_[j] = str;
        if(index_[j].m() != 1) 
            {
            str *= index_[j].m();
            ++rn_;
            }
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
    return IndexT();
    }

//
// Given IndexSetT<IndexT> iset and IndexT I,
// return int j such that iset[j] == I.
// If not found, returns -1
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
    return -1;
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
    return IndexT();
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
    for(const auto& J : is)
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
