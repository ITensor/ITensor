//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEXSET_H
#define __ITENSOR_INDEXSET_H
#include <algorithm>
#include "itensor/util/readwrite.h"
#include "itensor/util/safe_ptr.h"
#include "itensor/tensor/range.h"
#include "itensor/tensor/types.h"
#include "itensor/tensor/permutation.h"
#include "itensor/indexset_iter.h"
#include "itensor/index.h"


namespace itensor {

template<class IndexT>
class IndexSetT;

class IQIndex;

//
// IndexSetT
//
// Aliases:
using IndexSet = IndexSetT<Index>;
using IQIndexSet = IndexSetT<IQIndex>;

using IndexSetBuilder = RangeBuilderT<Index>;
using IQIndexSetBuilder = RangeBuilderT<IQIndex>;

//
// When constructed from a collection of indices,
// (as an explicit set of arguments or via
// a container) puts the indices with m>1 to the
// front and those with m==1 at the back but otherwise
// keeps the indices in the order given.
//

template <class index_type_>
class IndexSetT
    {
    public:
    using index_type = index_type_;
    using range_type = RangeT<index_type>;
    using size_type = typename range_type::size_type;
    using storage_type = typename range_type::storage_type;
    using value_type = index_type;
    using iterator = IndexSetIter<index_type>;
    using const_iterator = IndexSetIter<const index_type>;
    using indexval_type = typename index_type::indexval_type;
    private:
    range_type range_;
    public:

    IndexSetT() { }

    // construct from 1 or more indices
    template <typename... Inds>
    explicit
    IndexSetT(index_type const& i1, 
              Inds&&... rest)
      : range_(i1,std::forward<Inds>(rest)...)
        { }

    explicit
    IndexSetT(std::vector<index_type> const& ii) : range_(ii) { }

    template<size_t N>
    explicit
    IndexSetT(std::array<index_type,N> const& ii) : range_(ii) { }

    IndexSetT(std::initializer_list<index_type> ii) : range_(ii) { }

    explicit
    IndexSetT(RangeBuilderT<index_type> & builder) 
      : range_(range_type(builder)) 
        { }

    IndexSetT&
    operator=(storage_type&& store)
        {
        range_ = std::move(store);
        return *this;
        }
    
    explicit operator bool() const { return !range_.empty(); }

    long
    extent(size_type i) const { return range_.extent(i); }

    size_type
    stride(size_type i) const { return range_.stride(i); }

    long
    r() const { return range_.r(); }
    
    size_t
    size() const { return range_.r(); }

    bool
    empty() const { return range_.empty(); }

    // 0-indexed access
    index_type const&
    operator[](size_type i) const 
        { 
#ifdef DEBUG
        if(i >= size()) Error("IndexSetT[i] arg out of range");
#endif
        return range_[i].ext; 
        }

    // 1-indexed access
    index_type const&
    index(size_type I) const 
        { 
#ifdef DEBUG
        if(I < 1 || I > size()) Error("IndexSetT.index(i) arg out of range");
#endif
        return range_[I-1].ext; 
        }

    index_type const&
    front() const { return range_.front().ext; }

    index_type const&
    back() const { return range_.back().ext; }

    iterator
    begin() { return iterator{range_.data()}; }

    iterator
    end() { return iterator{range_.data()+range_.size()}; }

    const_iterator
    begin() const { return const_iterator{range_.data()}; }

    const_iterator
    end() const { return const_iterator{range_.data()+range_.size()}; }

    const_iterator
    cbegin() const { return begin(); }

    const_iterator
    cend() const { return end(); }

    void
    swap(IndexSetT& other) { range_.swap(other.range_); }

    void
    clear() { range_.clear(); }

    void
    dag() { for(auto& J : *this) J.dag(); }

    void
    read(std::istream& s);

    void
    write(std::ostream& s) const;

    void
    computeStrides() { range_.computeStrides(); }
    };

//
// IndexSetT Primelevel Methods
//

template<typename IndexT, typename... Types>
void 
prime(IndexSetT<IndexT>& is, 
      IndexType type,
      int inc);

template<typename IndexT, typename... Types>
void 
prime(IndexSetT<IndexT>& is, 
      IndexType type1,
      Types&&... rest);

template<typename IndexT>
void 
prime(IndexSetT<IndexT>& is, int inc = 1) { prime(is,All,inc); }

template<typename IndexT, typename... IVals>
void 
prime(IndexSetT<IndexT>& is,
      const typename IndexT::indexval_type& iv1,
      IVals&&... ivs);

//
//Given a list of indices and an increment (an int)
//as the optional last argument (default is inc=1)
//increment all indices NOT listed in the arguments
//by the amout inc.
//
//For example, primeExcept(is,I1,I3,I4,I7,2);
//will increment all prime levels by 2 except for
//those of I1,I3,I4, and I7.
//
template<typename IndexT, typename... Inds>
void 
primeExcept(IndexSetT<IndexT>& is, 
            const IndexT& I1, 
            Inds&&... inds);

template<typename IndexT, typename... ITs>
void 
primeExcept(IndexSetT<IndexT>& is, 
            IndexType it,
            ITs&&... etc);

template<typename IndexT>
void 
noprime(IndexSetT<IndexT>& is, IndexType type = All);

template<typename IndexT, typename... ITs>
void 
noprime(IndexSetT<IndexT>& is,
        IndexType it1,
        IndexType it2,
        ITs&&... rest);

template<typename IndexT, typename... Inds>
void 
noprime(IndexSetT<IndexT>& is, 
        const IndexT& I1, 
        Inds&&... inds);

template<typename IndexT>
void 
mapprime(IndexSetT<IndexT>& is, int plevold, int plevnew, IndexType type = All);

//
// IndexSetT helper methods
//


template<class IndexT>
Arrow
dir(const IndexSetT<IndexT>& is, const IndexT& I);


template <class IndexT>
const IndexT&
finddir(const IndexSetT<IndexT>& iset, Arrow dir);

//
// Given IndexSetT<IndexT> iset and IndexT I,
// return int j such that iset[j] == I.
// If not found, returns -1
//
template <class IndexT>
long
findindex(const IndexSetT<IndexT>& iset, 
          const IndexT& I);

template <class IndexT>
const IndexT&
findtype(const IndexSetT<IndexT>& iset, IndexType t);

//
// Compute the permutation P taking an IndexSetT iset
// to oset (of type IndexSetT or array<IndexT,NMAX>)
//
template <class IndexT>
void
getperm(const IndexSetT<IndexT>& iset, 
        const typename IndexSetT<IndexT>::storage& oset, 
        Permutation& P);

template <class IndexT>
bool
hasindex(const IndexSetT<IndexT>& iset, 
         const IndexT& I);

template <class IndexT>
bool
hastype(const IndexSetT<IndexT>& iset, 
        IndexType t);

template <class IndexT>
long
minM(const IndexSetT<IndexT>& iset);

template <class IndexT>
long
maxM(const IndexSetT<IndexT>& iset);

template<class IndexT>
void
contractIS(IndexSetT<IndexT> const& Lis,
           IndexSetT<IndexT> const& Ris,
           IndexSetT<IndexT> & Nis,
           bool sortResult = false);

template<class IndexT, class ContainerT>
void
contractIS(IndexSetT<IndexT> const& Lis,
           ContainerT const& Lind,
           IndexSetT<IndexT> const& Ris,
           ContainerT const& Rind,
           IndexSetT<IndexT> & Nis,
           ContainerT & Nind,
           bool sortResult = false);

template <class IndexT>
std::ostream&
operator<<(std::ostream& s, IndexSetT<IndexT> const& is);

} //namespace itensor

#include "itensor/indexset.ih"

#endif
