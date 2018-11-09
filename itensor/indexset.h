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

//
// IndexSet Primelevel Methods
//

// increment primelevel of all
// indices by an amount "inc"
void 
prime(IndexSet& is, 
      int inc = 1);

template<typename... VArgs>
void
primeLevel(IndexSet& is,
           VArgs&&... vargs);

//
// Increment primelevels of indices in the
// set by matching them against a list
// of other objects, including:
// * Index (or IQIndex) objects
// * IndexType objects
// * IndexVal (or IQIndexVal) objects
// The last argument can optionally
// be an integer "inc" telling how
// much to increment by.
//
template<typename... VArgs>
void 
prime(IndexSet& is, 
      VArgs&&... vargs);

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
template<typename... Inds>
void 
primeExcept(IndexSet& is, 
            Index const& I1, 
            Inds&&... inds);

template<typename... ITs>
void 
primeExcept(IndexSet& is, 
            IndexType it,
            ITs&&... etc);

void 
noprime(IndexSet& is, IndexType type = All);

template<typename... ITs>
void 
noprime(IndexSet& is,
        IndexType it1,
        IndexType it2,
        ITs&&... rest);

template<typename... Inds>
void 
noprime(IndexSet& is, 
        Index const& I1, 
        Inds&&... inds);

// This version of mapprime takes
// any number of triples: I,p1,p2
// where I is an index or an IndexType,
// p1 is the inital primelevel and
// p2 if the final primlevel
// For example,
// prime(is,j,0,2,Site,1,0);
// would change an Index "j" with
// primelevel 0 to have primelevel 2
// and any indices with type Site
// and primelevel 1 to have primelevel 0
// If two or more mappings match for a particular
// index, only the first mapping is used.
template<typename... VArgs>
void 
mapprime(IndexSet& is, 
         VArgs&&... vargs);

void 
mapprime(IndexSet& is, 
         int plevold, 
         int plevnew, 
         IndexType type = All);

//Replace all indices of type t by 'similar' indices 
//with same properties but which don't compare equal 
//to the indices they replace (using sim(IndexT) function)
void 
sim(IndexSet & is, 
    IndexType t);

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
findindex(IndexSet const& iset, 
          Index const& I);

Index
findtype(IndexSet const& iset, 
         IndexType t);

Index
finddir(IndexSet const& iset, Arrow dir);

Arrow
dir(IndexSet const& is, Index const& I);


bool
hasindex(IndexSet const& iset, 
         Index const& I);

bool
hastype(IndexSet const& iset, 
        IndexType t);

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
