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

    explicit
    IndexSetT(std::vector<index_type> const& ii) : parent(ii) { }

    template<size_t N>
    explicit
    IndexSetT(std::array<index_type,N> const& ii) : parent(ii) { }

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
        if(i >= parent::size()) Error("IndexSetT[i] arg out of range");
#endif
        return parent::index(i);
        }

    // 1-indexed access
    index_type &
    index(size_type I)
        { 
#ifdef DEBUG
        if(I < 1 || I > parent::size()) Error("IndexSetT.index(i) arg out of range");
#endif
        return operator[](I-1);
        }

    // 0-indexed access
    index_type const&
    operator[](size_type i) const
        { 
#ifdef DEBUG
        if(i >= parent::size()) Error("IndexSetT[i] arg out of range");
#endif
        return parent::index(i);
        }

    // 1-indexed access
    index_type const&
    index(size_type I) const
        { 
#ifdef DEBUG
        if(I < 1 || I > parent::size()) Error("IndexSetT.index(i) arg out of range");
#endif
        return operator[](I-1);
        }


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

    parent const&
    range() const { return *this; }

    void
    swap(IndexSetT& other) { parent::swap(other); }

    void
    dag() { for(auto& J : *this) J.dag(); }

    void
    read(std::istream& s);

    void
    write(std::ostream& s) const;
    };

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
IndexT const&
finddir(IndexSetT<IndexT> const& iset, Arrow dir);

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
typename IndexSetIter<T>::difference_type 
operator-(const IndexSetIter<T>& x, const IndexSetIter<T>& y) 
    { 
    return x.data() - y.data();
    } 

template <typename T>
IndexSetIter<T> 
operator+(const IndexSetIter<T>& x, typename IndexSetIter<T>::difference_type d) 
    { 
    return x += d;
    } 

template <typename T>
IndexSetIter<T> 
operator+(typename IndexSetIter<T>::difference_type d, const IndexSetIter<T>& x) 
    { 
    return x += d;
    } 

} //namespace itensor

#include "itensor/indexset.ih"

#endif
