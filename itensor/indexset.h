//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEXSET_H
#define __ITENSOR_INDEXSET_H
#include <algorithm>
#include "itensor/util/safe_ptr.h"
#include "itensor/index.h"
#include "itensor/iqindex.h"
#include "itensor/tensor/contract.h"
#include "itensor/tensor/range.h"
#include "itensor/tensor/types.h"
#include "itensor/tensor/permutation.h"

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

using IndexSetBuilder = RangeBuilderT<IndexSet>;
using IQIndexSetBuilder = RangeBuilderT<IQIndexSet>;

//
// When constructed from a collection of indices,
// (as an explicit set of arguments or via
// a container) puts the indices with m>1 to the
// front and those with m==1 at the back but otherwise
// keeps the indices in the order given.
//

template<typename index_type_> 
class IndexSetIter;

template <class index_type_>
class IndexSetT : public RangeT<index_type_>
    {
    public:
    using index_type = index_type_;
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

    IndexSetT() { }

    // construct from 1 or more indices
    template <typename... Inds>
    explicit
    IndexSetT(index_type const& i1, 
              Inds&&... rest)
      : parent(i1,std::forward<Inds>(rest)...)
        { }

    template<typename IndxContainer>
    explicit
    IndexSetT(IndxContainer && ii) 
      : parent(std::forward<IndxContainer>(ii)) { }

    IndexSetT(std::initializer_list<index_type> ii) : parent(ii) { }

    explicit
    IndexSetT(storage_type && store) 
      : parent(std::move(store)) 
        { }

    IndexSetT&
    operator=(storage_type&& store)
        {
        parent::operator=(std::move(store));
        return *this;
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
        if(i >= parent::size()) throw ITError("IndexSetT[i] arg out of range");
#endif
        return parent::index(i);
        }

    // 1-indexed access
    index_type &
    index(size_type I)
        { 
#ifdef DEBUG
        if(I < 1 || I > parent::size()) throw ITError("IndexSetT.index(i) arg out of range");
#endif
        return operator[](I-1);
        }

    // 0-indexed access
    index_type const&
    operator[](size_type i) const
        { 
#ifdef DEBUG
        if(i >= parent::size()) throw ITError("IndexSetT[i] arg out of range");
#endif
        return parent::index(i);
        }

    // 1-indexed access
    index_type const&
    index(size_type I) const
        { 
#ifdef DEBUG
        if(I < 1 || I > parent::size()) throw ITError("IndexSetT.index(i) arg out of range");
#endif
        return operator[](I-1);
        }

    parent const&
    range() const { return *this; }

    void
    dag() { for(auto& J : *this) J.dag(); }

    void
    swap(IndexSetT & other) { parent::swap(other); }

    index_type const&
    front() const { return parent::front().ind; }

    index_type const&
    back() const { return parent::back().ind; }

    iterator
    begin() { return iterator{*this}; }

    iterator
    end() { return iterator::makeEnd(*this); }

    const_iterator
    begin() const { return const_iterator{*this}; }

    const_iterator
    end() const { return const_iterator::makeEnd(*this); }

    const_iterator
    cbegin() const { return begin(); }

    const_iterator
    cend() const { return end(); }

    };

template<typename index_type>
void
read(std::istream& s, IndexSetT<index_type> & is);

template<typename index_type>
void
write(std::ostream& s, IndexSetT<index_type> const& is);

template<typename index_type>
auto
rangeBegin(IndexSetT<index_type> const& is) -> decltype(is.range().begin())
    {
    return is.range().begin();
    }

template<typename index_type>
auto
rangeEnd(IndexSetT<index_type> const& is) -> decltype(is.range().end())
    {
    return is.range().end();
    }

//
// IndexSetT Primelevel Methods
//

// increment primelevel of all
// indices by an amount "inc"
template<typename IndexT>
void 
prime(IndexSetT<IndexT>& is, 
      int inc = 1);

template<typename IndexT, typename... VArgs>
void
primeLevel(IndexSetT<IndexT>& is,
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
template<typename IndexT, typename... VArgs>
void 
prime(IndexSetT<IndexT>& is, 
      VArgs&&... vargs);

//// Increment primelevels of the indices
//// specified by 1, or an optional amount "inc"
//// For example, to prime indices I and J by 2,
//// prime(is,I,J,2);
//template<typename IndexT, typename... Inds>
//void 
//prime(IndexSetT<IndexT>& is, 
//      IndexT const& I1, 
//      Inds&&... rest);

//// increment primelevel of all indices of
//// type "type" by an amount "inc"
//template<typename IndexT, typename... Types>
//void 
//prime(IndexSetT<IndexT>& is, 
//      IndexType type,
//      int inc = 1);

//// same as above but for multiple types
//// optionally, last argument can be 
//// an increment amount
//template<typename IndexT, typename... Types>
//void 
//prime(IndexSetT<IndexT>& is, 
//      IndexType type1,
//      Types&&... rest);

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
            IndexT const& I1, 
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
        IndexT const& I1, 
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
template<typename IndexT, typename... VArgs>
void 
mapprime(IndexSetT<IndexT>& is, 
         VArgs&&... vargs);

template<typename IndexT>
void 
mapprime(IndexSetT<IndexT>& is, 
         int plevold, 
         int plevnew, 
         IndexType type = All);

//Replace all indices of type t by 'similar' indices 
//with same properties but which don't compare equal 
//to the indices they replace (using sim(IndexT) function)
template<typename IndexT>
void 
sim(IndexSetT<IndexT> & is, 
    IndexType t);

//Replace index I with a 'similar' index having same properties
//but which does not compare equal to it (using sim(I) function)
template<typename IndexT>
void 
sim(IndexSetT<IndexT> & is, 
    IndexT const& I);

//
// IndexSetT helper methods
//


//
// Given IndexSetT<IndexT> iset and IndexT I,
// return int j such that iset[j] == I.
// If not found, returns -1
//
template <class IndexT>
long
findindex(IndexSetT<IndexT> const& iset, 
          IndexT const& I);

template <class IndexT>
IndexT const&
findtype(IndexSetT<IndexT> const& iset, 
         IndexType t);

template <class IndexT>
IndexT const&
finddir(IndexSetT<IndexT> const& iset, Arrow dir);

template<class IndexT>
Arrow
dir(IndexSetT<IndexT> const& is, IndexT const& I);


////
//// Compute the permutation P taking an IndexSetT iset
//// to oset (of type IndexSetT or array<IndexT,NMAX>)
////
//template <class IndexT>
//void
//getperm(const IndexSetT<IndexT>& iset, 
//        const typename IndexSetT<IndexT>::storage& oset, 
//        Permutation& P);

template <class IndexT>
bool
hasindex(IndexSetT<IndexT> const& iset, 
         IndexT const& I);

template <class IndexT>
bool
hastype(IndexSetT<IndexT> const& iset, 
        IndexType t);

template <class IndexT>
long
minM(IndexSetT<IndexT> const& iset);

template <class IndexT>
long
maxM(IndexSetT<IndexT> const& iset);

template<class IndexT>
void
contractIS(IndexSetT<IndexT> const& Lis,
           IndexSetT<IndexT> const& Ris,
           IndexSetT<IndexT> & Nis,
           bool sortResult = false);

template<class IndexT, class LabelT>
void
contractIS(IndexSetT<IndexT> const& Lis,
           LabelT const& Lind,
           IndexSetT<IndexT> const& Ris,
           LabelT const& Rind,
           IndexSetT<IndexT> & Nis,
           LabelT & Nind,
           bool sortResult = false);

template<class IndexT, class LabelT>
void
ncprod(IndexSetT<IndexT> const& Lis,
       LabelT const& Lind,
       IndexSetT<IndexT> const& Ris,
       LabelT const& Rind,
       IndexSetT<IndexT> & Nis,
       LabelT & Nind);

template <class IndexT>
std::ostream&
operator<<(std::ostream& s, IndexSetT<IndexT> const& is);

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
                                             const IndexSetT<index_type>,
                                             IndexSetT<index_type>>;
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

template<typename IndexT>
struct IndexValIter
    {
    using indexval_type = typename IndexT::indexval_type;

    IndexSetT<IndexT> const& is;
    detail::GCounter count;
    bool done = false;
    IndexValIter(IndexSetT<IndexT> const& is_) 
      : is(is_),
        count(is_.size())
        { 
        for(auto n : range(is.size()))
            {
            count.setRange(n,0,is[n].m()-1);
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
        //    eit.count.setRange(n,is[n].m()-1,is[n].m()-1);
        //    }
        return eit;
        }

    bool
    operator!=(IndexValIter const& other) { return done != other.done; }

    std::vector<indexval_type>
    operator*() 
        { 
        auto res = std::vector<indexval_type>(is.size());
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

template<class I>
detail::IndexValIter<I>
iterInds(IndexSetT<I> const& is)
    {
    return detail::IndexValIter<I>(is);
    }

} //namespace itensor

#include "itensor/indexset_impl.h"

#endif
