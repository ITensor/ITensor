//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_INDEXSET_H
#define __ITENSOR_INDEXSET_H
#include <algorithm>
#include "itensor/util/safe_ptr.h"
#include "itensor/index.h"
#include "itensor/tensor/contract.h"
#include "itensor/tensor/range.h"
#include "itensor/tensor/types.h"
#include "itensor/tensor/permutation.h"

namespace itensor {


class IndexSet;

using IndexSetBuilder = RangeBuilderT<IndexSet>;

template<typename I>
class IndexSetIter;

//
// IndexSet
//
// When constructed from a collection of indices,
// (as an explicit set of arguments or via
// a container) puts the indices with m>1 to the
// front and those with m==1 at the back but otherwise
// keeps the indices in the order given.
//

void
checkQNConsistent(IndexSet const&);

class IndexSet : public RangeT<Index>
    {
    public:
    using extent_type = index_type;
    using range_type = RangeT<Index>;
    using parent = RangeT<Index>;
    using size_type = typename range_type::size_type;
    using storage_type = typename range_type::storage_type;
    using value_type = Index;
    using iterator = IndexSetIter<Index>;
    using const_iterator = IndexSetIter<const Index>;

    public:

    IndexSet() { }

    // construct from 1 or more indices
    template <typename... Inds>
    explicit
    IndexSet(Index const& i1, 
             Inds&&... inds)
      : parent(i1,std::forward<Inds>(inds)...)
        { 
        checkQNConsistent(*this);
        }

    IndexSet(std::initializer_list<Index> const& ii)
      : parent(ii) 
        { 
        checkQNConsistent(*this);
        }

    IndexSet(std::vector<Index> const& ii)
      : parent(ii) 
        { 
        checkQNConsistent(*this);
        }

    template<size_t N>
    IndexSet(std::array<Index,N> const& ii)
      : parent(ii)
        {
        checkQNConsistent(*this);
        }

    template<typename IndxContainer>
    explicit
    IndexSet(IndxContainer && ii) 
      : parent(std::forward<IndxContainer>(ii)) 
        { 
        checkQNConsistent(*this);
        }

    explicit
    IndexSet(storage_type && store) 
      : parent(std::move(store)) 
        { 
        checkQNConsistent(*this);
        }

    // construct from 2 IndexSets
    IndexSet(IndexSet const& is1,
             IndexSet const& is2)
        {
        auto N1 = is1.order();
        auto N2 = is2.order();
        auto N = N1+N2;
        auto inds = IndexSetBuilder(N);
        for( auto n1 : range1(N1) )
          inds.nextIndex(std::move(is1(n1)));
        for( auto n2 : range1(N2) )
          inds.nextIndex(std::move(is2(n2)));
        *this = inds.build();
        checkQNConsistent(*this);
        }

    // construct from an Index and IndexSet
    IndexSet(Index const& i,
             IndexSet const& is)
        {
        auto N = is.order();
        auto inds = IndexSetBuilder(N+1);
        inds.nextIndex(std::move(i));
        for( auto n : range1(N) )
          inds.nextIndex(std::move(is(n)));
        *this = inds.build();
        checkQNConsistent(*this);
        }

    // construct from an Index and IndexSet
    IndexSet(IndexSet const& is,
             Index const& i)
        {
        auto N = is.order();
        auto inds = IndexSetBuilder(N+1);
        for( auto n : range1(N) )
          inds.nextIndex(std::move(is(n)));
        inds.nextIndex(std::move(i));
        *this = inds.build();
        checkQNConsistent(*this);
        }

    explicit operator bool() const { return !parent::empty(); }

    long
    extent(size_type i) const { return parent::extent(i); }

    size_type
    stride(size_type i) const { return parent::stride(i); }

    // Get the number of indices
    long
    order() const { return parent::order(); }
 
    // Get the number of indices (alternative)
    long
    length() const { return this->order(); }
    
    // 0-indexed access
    Index &
    operator[](size_type i)
        { 
#ifdef DEBUG
        if(i >= parent::size()) throw ITError("IndexSet[i] arg out of range");
#endif
        return parent::index(i);
        }

    // 1-indexed access
    Index &
    operator()(size_type I)
        { 
#ifdef DEBUG
        if(I < 1 || I > parent::size()) throw ITError("IndexSet(i) arg out of range");
#endif
        return operator[](I-1);
        }

    // Deprecated
    Index &
    index(size_type I) { return operator()(I); }

    // 0-indexed access
    Index const&
    operator[](size_type i) const
        { 
#ifdef DEBUG
        if(i >= parent::size()) throw ITError("IndexSet[i] arg out of range");
#endif
        return parent::index(i);
        }

    // 1-indexed access
    Index const&
    operator()(size_type I) const
        { 
#ifdef DEBUG
        if(I < 1 || I > parent::size()) throw ITError("IndexSet(i) arg out of range");
#endif
        return operator[](I-1);
        }

    // Deprecated
    Index const&
    index(size_type I) const { return operator()(I); }

    parent const&
    range() const { return *this; }

    IndexSet&
    dag();

    void
    swap(IndexSet & other) { parent::swap(other); }

    Index const&
    front() const { return parent::front().ind; }

    Index const&
    back() const { return parent::back().ind; }

    iterator
    begin();

    iterator
    end();

    const_iterator
    begin() const;

    const_iterator
    end() const;

    const_iterator
    cbegin() const;

    const_iterator
    cend() const;

    //
    // Tag methods
    //

    IndexSet&
    setTags(TagSet const& tsnew);

    IndexSet&
    setTags(TagSet const& tsnew, 
            IndexSet const& ismatch);

    template<typename... VarArgs>
    IndexSet&
    setTags(TagSet const& tsnew,
            Index const& imatch1,
            VarArgs&&... vargs)
      {
      setTags(tsnew,IndexSet(imatch1,std::forward<VarArgs>(vargs)...));
      return *this;
      }

    IndexSet&
    setTags(TagSet const& tsnew, 
            TagSet const& tsmatch);

    IndexSet&
    noTags();

    IndexSet&
    noTags(IndexSet const& ismatch);

    template<typename... VarArgs>
    IndexSet&
    noTags(Index const& imatch1,
           VarArgs&&... vargs)
      {
      noTags(IndexSet(imatch1,std::forward<VarArgs>(vargs)...));
      return *this;
      }

    IndexSet&
    noTags(TagSet const& tsmatch);

    IndexSet&
    addTags(TagSet const& tsadd);

    IndexSet&
    addTags(TagSet const& tsadd, 
            IndexSet const& ismatch);

    template<typename... VarArgs>
    IndexSet&
    addTags(TagSet const& tsadd,
            Index const& imatch1,
            VarArgs&&... vargs)
      {
      addTags(tsadd,IndexSet(imatch1,std::forward<VarArgs>(vargs)...));
      return *this;
      }

    IndexSet&
    addTags(TagSet const& tsadd,
            TagSet const& tsmatch);

    IndexSet&
    removeTags(TagSet const& tsremove);

    IndexSet&
    removeTags(TagSet const& tsremove, 
               IndexSet const& ismatch);

    template<typename... VarArgs>
    IndexSet&
    removeTags(TagSet const& tsremove,
               Index const& imatch1,
               VarArgs&&... vargs)
      {
      removeTags(tsremove,IndexSet(imatch1,std::forward<VarArgs>(vargs)...));
      return *this;
      }

    IndexSet&
    removeTags(TagSet const& tsremove, 
               TagSet const& tsmatch);

    IndexSet&
    replaceTags(TagSet const& tsold, 
                TagSet const& tsnew);

    IndexSet&
    replaceTags(TagSet const& tsold, 
                TagSet const& tsnew, 
                IndexSet const& ismatch);

    template<typename... VarArgs>
    IndexSet&
    replaceTags(TagSet const& tsold,
                TagSet const& tsnew,
                Index const& imatch1,
                VarArgs&&... vargs)
      {
      replaceTags(tsold,tsnew,IndexSet(imatch1,std::forward<VarArgs>(vargs)...));
      return *this;
      }

    IndexSet&
    replaceTags(TagSet const& tsold, 
                TagSet const& tsnew,
                TagSet const& tsmatch);

    template<typename... VarArgs>
    IndexSet&
    swapTags(TagSet const& ts1,
             TagSet const& ts2,
             VarArgs&&... vargs);

    //
    // Integer tag convenience functions
    //

    //
    // Set the integer tag of indices to plnew
    //

    IndexSet&
    setPrime(int plnew);

    IndexSet&
    setPrime(int plnew,
             IndexSet const& ismatch);

    template<typename... VarArgs>
    IndexSet&
    setPrime(int plnew,
             Index const& imatch1,
             VarArgs&&... vargs)
        {
        setPrime(plnew,IndexSet(imatch1,std::forward<VarArgs>(vargs)...));
        return *this;
        }

    IndexSet&
    setPrime(int plnew,
             TagSet const& tsmatch);

    IndexSet&
    mapPrime(int plold, int plnew);

    IndexSet&
    mapPrime(int plold, int plnew,
             IndexSet const& ismatch);

    template<typename... VarArgs>
    IndexSet&
    mapPrime(int plold, int plnew,
             Index const& imatch1,
             VarArgs&&... vargs)
        {
        mapPrime(plold,plnew,IndexSet(imatch1,std::forward<VarArgs>(vargs)...));
        return *this;
        }

    IndexSet&
    mapPrime(int plold, int plnew,
             TagSet const& tsmatch);

    template<typename... VarArgs>
    IndexSet&
    swapPrime(int pl1,
              int pl2,
              VarArgs&&... vargs);

    template<typename... VarArgs>
    IndexSet&
    noPrime(VarArgs&&... vargs)
        {
        setPrime(0,std::forward<VarArgs>(vargs)...);
        return *this;
        }

    //
    // Increase the integer tag of indices by plinc
    //

    IndexSet&
    prime(int plinc);

    IndexSet&
    prime()
      {
      prime(1);
      return *this;
      }

    IndexSet&
    prime(int plinc,
          IndexSet const& ismatch);

    IndexSet&
    prime(IndexSet const& ismatch)
      {
      prime(1,ismatch);
      return *this;
      }

    template<typename... VarArgs>
    IndexSet&
    prime(int plinc,
          Index const& imatch1,
          VarArgs&&... vargs)
      {
      prime(plinc,IndexSet(imatch1,std::forward<VarArgs>(vargs)...));
      return *this;
      }

    template<typename... VarArgs>
    IndexSet&
    prime(Index const& imatch1,
          VarArgs&&... vargs)
      {
      prime(IndexSet(imatch1,std::forward<VarArgs>(vargs)...));
      return *this;
      }

    IndexSet&
    prime(int plinc,
          TagSet const& tsmatch);

    IndexSet&
    prime(TagSet const& tsmatch)
      {
      prime(1,tsmatch);
      return *this;
      }
 
    // Remove QNs from all indices in the IndexSet
    IndexSet&
    removeQNs();

    //
    // Deprecated
    //

    long
    r() const { return this->order(); }
    
    void
    prime(Index const& imatch,
          int plinc)
        {
        Error("Error: .prime(Index,int) is no longer supported, use .prime(int,Index) instead.");
        }

    };

void
read(std::istream& s, IndexSet & is);

void
write(std::ostream& s, IndexSet const& is);

auto inline
rangeBegin(IndexSet const& is) -> decltype(is.range().begin())
    {
    return is.range().begin();
    }

auto inline
rangeEnd(IndexSet const& is) -> decltype(is.range().end())
    {
    return is.range().end();
    }

long
order(IndexSet const& is);

long
length(IndexSet const& is);

IndexSet
dag(IndexSet is);

template<typename... VarArgs>
IndexSet
setTags(IndexSet A,
        VarArgs&&... vargs)
    {
    A.setTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
IndexSet
noTags(IndexSet A,
       VarArgs&&... vargs)
    {
    A.noTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
IndexSet
addTags(IndexSet A,
        VarArgs&&... vargs)
    {
    A.addTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
IndexSet
removeTags(IndexSet A,
           VarArgs&&... vargs)
    {
    A.removeTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
IndexSet
replaceTags(IndexSet A,
            VarArgs&&... vargs)
    {
    A.replaceTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
IndexSet
swapTags(IndexSet A,
         VarArgs&&... vargs)
    {
    A.swapTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
IndexSet
prime(IndexSet A,
      VarArgs&&... vargs)
    {
    A.prime(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
IndexSet
setPrime(IndexSet A,
         VarArgs&&... vargs)
    {
    A.setPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
IndexSet
mapPrime(IndexSet A,
         VarArgs&&... vargs)
    {
    A.mapPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
IndexSet
swapPrime(IndexSet A,
          VarArgs&&... vargs)
    {
    A.swapPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
IndexSet
noPrime(IndexSet A,
        VarArgs&&... vargs)
    {
    A.noPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

//Replace all indices with 'similar' indices 
//with the same properties but which don't compare equal 
//to the indices they replace (using sim(Index) function)
IndexSet
sim(IndexSet is);
IndexSet
sim(IndexSet is, 
    IndexSet const& ismatch);
IndexSet
sim(IndexSet is, 
    TagSet const& tsmatch);


//
// IndexSet helper methods
//


//
// Given IndexSet iset and Index I,
// return int j such that iset[j] == I.
// If not found, returns -1
//
int
indexPosition(IndexSet const& is, 
              Index const& imatch);

std::vector<int>
indexPositions(IndexSet const& is,
               IndexSet const& ismatch);

Arrow
dir(IndexSet const& is, Index const& I);

// Return true if the Index `imatch` is in
// `is`
bool
hasIndex(IndexSet const& is, 
         Index const& imatch);

// Return true if all indices of `ismatch`
// are in `is`
bool
hasInds(IndexSet const& is,
        IndexSet const& ismatch);

template<typename... Inds>
bool
hasInds(IndexSet const& is,
        Index const& i1, Inds&&... inds)
  {
  return hasInds(is,IndexSet(i1,std::forward<Inds>(inds)...));
  }

// Return true if IndexSet `is1` and `is2` have 
// the same indices
bool
hasSameInds(IndexSet const& is1,
            IndexSet const& is2);

// IndexSets are equal if they are the
// same size and contain equal indices
// in equal ordering (i.e. equals({i,j},{i,j}) -> true)
// but equals({i,j},{j,i}) -> false).
// For set equality, you can use hasSameInds(is1,is2).
bool
equals(IndexSet const& is1,
       IndexSet const& is2);

long
minDim(IndexSet const& iset);

long
maxDim(IndexSet const& iset);

void
contractIS(IndexSet const& Lis,
           IndexSet const& Ris,
           IndexSet & Nis,
           bool sortResult = false);

template<class LabelT>
void
contractIS(IndexSet const& Lis,
           LabelT const& Lind,
           IndexSet const& Ris,
           LabelT const& Rind,
           IndexSet & Nis,
           LabelT & Nind,
           bool sortResult = false);

template<class LabelT>
void
contractISReplaceIndex(IndexSet const& Lis,
                       LabelT const& Lind,
                       IndexSet const& Ris,
                       LabelT const& Rind,
                       IndexSet & Nis);

template<class LabelT>
void
ncprod(IndexSet const& Lis,
       LabelT const& Lind,
       IndexSet const& Ris,
       LabelT const& Rind,
       IndexSet & Nis,
       LabelT & Nind);

std::ostream&
operator<<(std::ostream& s, IndexSet const& is);

template<typename index_type_>
class IndexSetIter
    { 
    public:
    using index_type = stdx::remove_const_t<index_type_>;
    using value_type = index_type;
    using reference = index_type_&;
    using difference_type = std::ptrdiff_t;
    using pointer = index_type_*;
    using iterator_category = std::random_access_iterator_tag;
    using indexset_type = stdx::conditional_t<std::is_const<index_type_>::value,
                                             const IndexSet,
                                             IndexSet>;
    using range_ptr = typename RangeT<index_type>::value_type*;
    using const_range_ptr = const typename RangeT<index_type>::value_type*;
    using data_ptr = stdx::conditional_t<std::is_const<index_type_>::value,
                                      const_range_ptr,
                                      range_ptr>;
    private:
    size_t off_ = 0;
    indexset_type* p_; 
    public: 

    IndexSetIter() : p_(nullptr) { }

    explicit
    IndexSetIter(indexset_type & is) : p_(&is) { }

    size_t
    offset() const { return off_; }

    IndexSetIter& 
    operator++() 
        { 
        ++off_; 
        return *this; 
        } 

    IndexSetIter 
    operator++(int) 
        { 
        auto tmp = *this; //save copy of this
        ++off_; 
        return tmp; 
        } 

    IndexSetIter& 
    operator+=(difference_type x) 
        { 
        off_ += x;
        return *this; 
        } 

    IndexSetIter& 
    operator--( ) 
        { 
        --off_;
        return *this; 
        } 

    IndexSetIter 
    operator--(int) 
        { 
        auto tmp = *this; //save copy of this
        --off_;
        return tmp; 
        } 

    IndexSetIter& 
    operator-=(difference_type x) 
        { 
        off_ -= x;
        return *this; 
        } 

    reference 
    operator[](difference_type n) { return p_->operator[](n); } 

    reference 
    operator*() { return p_->operator[](off_); }  

    pointer 
    operator->() { return &(p_->operator[](off_)); }

    IndexSetIter static
    makeEnd(indexset_type & is)
        {
        IndexSetIter end;
        end.p_ = &is;
        end.off_ = is.size();
        return end;
        }
    }; 

template <typename T>
bool 
operator==(const IndexSetIter<T>& x, const IndexSetIter<T>& y) 
    { 
    return x.offset() == y.offset(); 
    } 

template <typename T>
bool 
operator!=(const IndexSetIter<T>& x, const IndexSetIter<T>& y) 
    { 
    return x.offset() != y.offset(); 
    } 

template <typename T>
bool 
operator<(const IndexSetIter<T>& x, const IndexSetIter<T>& y) 
    { 
    return x.offset() < y.offset(); 
    } 

template <typename T>
IndexSetIter<T>
operator+(IndexSetIter<T> x, 
          typename IndexSetIter<T>::difference_type d) 
    { 
    return x += d;
    } 

template <typename T>
IndexSetIter<T>
operator+(typename IndexSetIter<T>::difference_type d, 
          IndexSetIter<T> x) 
    { 
    return x += d;
    } 


//
// IndexValIter - helper for iterInds
//

namespace detail {

struct IndexValIter
    {
    IndexSet const& is;
    detail::GCounter count;
    bool done = false;
    IndexValIter(IndexSet const& is_) 
      : is(is_),
        count(is_.size())
        { 
        for(auto n : range(is.size()))
            {
            count.setRange(n,0,is[n].dim()-1);
            }
        }

    IndexValIter
    begin() const { return *this; }

    IndexValIter
    end() const 
        { 
        auto eit = *this;
        eit.done = true;
        //for(auto n : range(is.size())) 
        //    {
        //    eit.count.setRange(n,is[n].dim()-1,is[n].dim()-1);
        //    }
        return eit;
        }

    bool
    operator!=(IndexValIter const& other) { return done != other.done; }

    std::vector<IndexVal>
    operator*() 
        { 
        auto res = std::vector<IndexVal>(is.size());
        for(auto n : range(is.size())) 
            {
            res.at(n) = is[n](1+count[n]);
            }
        return res;
        }

    IndexValIter&
    operator++()
        {
        ++count;
        done = !count.notDone();
        return *this;
        }
    };

} //namespace detail

detail::IndexValIter
iterInds(IndexSet const& is);

bool
hasQNs(IndexSet const& is);

QN
flux(std::vector<IndexVal> const& ivs);

void
checkIndexSet(IndexSet const& is);

void
checkIndexPositions(std::vector<int> const& is);

//
// IndexSet set operations
//

// Find the Indices containing tags in the specified TagSet 
IndexSet
findInds(IndexSet const& is,
         TagSet const& tsmatch);

// Convert an order one IndexSet into an Index
// Throws an error if more than one Index is in the IndexSet
// If no indices are found, returns a null Index
Index
findIndex(IndexSet const& is);

// Find the Index containing tags in the specified TagSet 
// Throws an error if more than one Index is found
// If no indices are found, returns a null Index
Index
findIndex(IndexSet const& is,
          TagSet const& tsmatch);

// Find the Indices not containing tags in the specified TagSet 
// Same as uniqueInds(is,findInds(is,tsmatch))
IndexSet
findIndsExcept(IndexSet const& is,
               TagSet const& tsmatch);

// Intersection of two IndexSets
IndexSet
commonInds(IndexSet const& is1,
           IndexSet const& is2);

// Union of two IndexSets (is1+is2)
// Preserves the ordering of the original
// IndexSets
IndexSet
unionInds(IndexSet const& is1,
          IndexSet const& is2);

IndexSet
unionInds(Index const& i,
          IndexSet const& is);

IndexSet
unionInds(IndexSet const& is,
          Index const& i);

IndexSet
unionInds(std::vector<IndexSet> const& is1);

// Difference of two IndexSets (is1-is2)
IndexSet
uniqueInds(IndexSet const& is1,
           IndexSet const& is2);

// Difference of IndexSet from a set
// of other IndexSets (is1-(is2+is3+...))
IndexSet
uniqueInds(IndexSet const& is1,
           std::vector<IndexSet> const& is2);

// Symmetric difference of two IndexSets
// (union(is1-is2,is2-is1))
IndexSet
noncommonInds(IndexSet const& is1,
              IndexSet const& is2);

#ifdef ITENSOR_USE_HDF5
void
h5_write(h5::group parent, std::string const& name, IndexSet const& is);
void
h5_read(h5::group parent, std::string const& name, IndexSet & is);
#endif

} //namespace itensor

#include "itensor/indexset_impl.h"

#endif
