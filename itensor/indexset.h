//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
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
    using index_type = Index;
    using extent_type = index_type;
    using range_type = RangeT<index_type>;
    using parent = RangeT<index_type>;
    using size_type = typename range_type::size_type;
    using storage_type = typename range_type::storage_type;
    using value_type = index_type;
    using iterator = IndexSetIter<index_type>;
    using const_iterator = IndexSetIter<const index_type>;
    using indexval_type = typename index_type::indexval_type;

    public:

    IndexSet() { }

    // construct from 1 or more indices
    template <typename... Inds>
    explicit
    IndexSet(index_type const& i1, 
             Inds&&... rest)
      : parent(i1,std::forward<Inds>(rest)...)
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

    IndexSet(std::initializer_list<index_type> ii) : parent(ii) 
        { 
        checkQNConsistent(*this);
        }

    explicit
    IndexSet(storage_type && store) 
      : parent(std::move(store)) 
        { 
        checkQNConsistent(*this);
        }

    explicit operator bool() const { return !parent::empty(); }

    long
    extent(size_type i) const { return parent::extent(i); }

    size_type
    stride(size_type i) const { return parent::stride(i); }

    long
    r() const { return parent::r(); }
    
    // 0-indexed access
    index_type &
    operator[](size_type i)
        { 
#ifdef DEBUG
        if(i >= parent::size()) throw ITError("IndexSet[i] arg out of range");
#endif
        return parent::index(i);
        }

    // 1-indexed access
    index_type &
    index(size_type I)
        { 
#ifdef DEBUG
        if(I < 1 || I > parent::size()) throw ITError("IndexSet.index(i) arg out of range");
#endif
        return operator[](I-1);
        }

    // 0-indexed access
    index_type const&
    operator[](size_type i) const
        { 
#ifdef DEBUG
        if(i >= parent::size()) throw ITError("IndexSet[i] arg out of range");
#endif
        return parent::index(i);
        }

    // 1-indexed access
    index_type const&
    index(size_type I) const
        { 
#ifdef DEBUG
        if(I < 1 || I > parent::size()) throw ITError("IndexSet.index(i) arg out of range");
#endif
        return operator[](I-1);
        }

    parent const&
    range() const { return *this; }

    void
    dag();

    void
    swap(IndexSet & other) { parent::swap(other); }

    index_type const&
    front() const { return parent::front().ind; }

    index_type const&
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
    // Prime methods
    //
    void
    setPrime(int plnew,
             TagSet const& tsmatch = TagSet(All));

    void
    setPrime(int plnew, Index const& imatch);

    template<typename... VarArgs>
    void
    setPrime(int plnew1, Index const& imatch1,
             int plnew2, Index const& imatch2,
             VarArgs&&... vargs);

    template<typename... VarArgs>
    void
    setPrime(int plnew,
             Index const& imatch1,
             Index const& imatch2,
             VarArgs&&... vargs);

    void
    noPrime(TagSet const& tsmatch = TagSet(All));

    void
    noPrime(Index const& imatch);

    //template<typename... VarArgs>
    //void
    //noPrime(Index const& imatch1,
    //        Index const& imatch2,
    //        VarArgs&&... vargs);

    void
    mapPrime(int plold, int plnew,
             TagSet const& tsmatch = TagSet(All));

    void
    prime(int plinc = 1,
          TagSet const& tsmatch = TagSet(All));

    void
    prime(TagSet const& tsmatch);

    void
    prime(int plinc, Index const& imatch);

    void
    prime(Index const& imatch);

    void
    swapPrime(int pl1, int pl2,
              TagSet const& tsmatch = TagSet());

    //
    // Tag methods
    //

    //TODO: add replaceTags(tsold,tsnew,index1,index2,...)?
    void
    replaceTags(TagSet const& tsold, 
                TagSet const& tsnew, 
                TagSet const& tsmatch = TagSet(All),
                int plmatch = -1);

    //TODO: add swapTags(ts1,ts2,index1,index2,...)?
    void
    swapTags(TagSet const& ts1, 
             TagSet const& ts2, 
             TagSet const& tsmatch = TagSet(All), 
             int plmatch = -1);

    //TODO: add setTags(tsnew1,imatch1,tsnew2,imatch2,...)?
    void
    setTags(TagSet const& tsnew, 
            TagSet const& tsmatch = TagSet(All), 
            int plmatch = -1);

    void
    addTags(TagSet const& tsadd, 
            TagSet const& tsmatch = TagSet(All), 
            int plmatch = -1);

    void
    addTags(TagSet const& tsadd, 
            int plmatch);

    //TODO: add addTags(tsadd1,imatch1,tsadd2,imatch2,...)?
    void
    addTags(TagSet const& tsadd, 
            Index const& imatch);

    //TODO: add removeTags(tsremove = TagSet(All),...)?
    void
    removeTags(TagSet const& tsremove, 
               TagSet const& tsmatch = TagSet(All), 
               int plmatch = -1);

    void
    removeTags(TagSet const& tsremove, 
               int plmatch);

    //TODO: add removeTags(tsremove1,imatch1,tsremove2,imatch2,...)?
    void
    removeTags(TagSet const& tsremove, 
               Index const& imatch);

    // Remove QNs from all indices in the IndexSet
    void
    removeQNs();

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

// Return the first index found with the specified tags
// This is useful if we know there is only one index
// that contains Tags ts, but don't know the other tags
Index
findIndex(IndexSet const& is,
          TagSet const& tsmatch, 
          int plmatch = -1);

// Find the Index with a certain TagSet and prime level
// If multiple indices or no index is found, throw an error
Index
findIndexExact(IndexSet const& is,
               TagSet const& tsmatch, 
               int plmatch = -1);
//
//
// IndexSet Primelevel Methods
//


//Replace all indices of type t by 'similar' indices 
//with same properties but which don't compare equal 
//to the indices they replace (using sim(Index) function)
void 
sim(IndexSet & is, 
    TagSet const& t);

//Replace index I with a 'similar' index having same properties
//but which does not compare equal to it (using sim(I) function)
void 
sim(IndexSet & is, 
    Index const& I);

//
// IndexSetT helper methods
//


//
// Given IndexSetT iset and Index I,
// return int j such that iset[j] == I.
// If not found, returns -1
//
long
findIndex(IndexSet const& iset, 
          Index const& I);

Index
finddir(IndexSet const& iset, Arrow dir);

Arrow
dir(IndexSet const& is, Index const& I);


bool
hasIndex(IndexSet const& iset, 
         Index const& I);

long
minM(IndexSet const& iset);

long
maxM(IndexSet const& iset);

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

bool
hasQNs(IndexSet const& is);

bool
hasQNs(std::vector<Index> const& inds);

} //namespace itensor

#include "itensor/indexset_impl.h"

#endif
