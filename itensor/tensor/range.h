//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_RANGE_H
#define __ITENSOR_RANGE_H

#include <array>
#include <vector>
#include "itensor/util/readwrite.h"
#include "itensor/util/autovector.h"
#include "itensor/util/vararray.h"
#include "itensor/util/infarray.h"

namespace itensor {

template<typename extent_type_>
class RangeT;

using Range = RangeT<size_t>;

//Storage type for RangeT
template<typename extent_type>
struct ExtStr
    {
    using size_type = size_t;
    extent_type ext = extent_type{}; //convertible to size_type
    size_type   str = 0; //stride

    ExtStr(extent_type x, size_type s) : ext(x), str(s) { }

    ExtStr() { }

    void
    write(std::ostream& s) const
        {
        itensor::write(s,ext);
        itensor::write(s,str);
        }
    void
    read(std::istream& s)
        {
        itensor::read(s,ext);
        itensor::read(s,str);
        }
    };

template<typename extent_type_>
class RangeT
    {
    public:
    using size_type = size_t;
    using extent_type = extent_type_;
    using value_type = ExtStr<extent_type>;
    using storage_type = std::vector<value_type>;
    private:
    storage_type store_;
    public:

    RangeT() { }

    template <typename... Inds>
    explicit
    RangeT(extent_type i0, 
           Inds&&... rest)
        {
        std::array<extent_type,1+sizeof...(rest)> inds = 
            {{i0,static_cast<extent_type>(rest)...}};
        init(inds);
        }

    explicit
    RangeT(const std::vector<extent_type>& v)
        {
        init(v);
        }

    template<size_t size>
    explicit
    RangeT(const std::array<extent_type,size>& a)
        {
        init(a);
        }

    RangeT(std::initializer_list<extent_type> ii) { init(ii); }

    //Most efficient way to construct range is by
    //providing a storage_type object with ext's
    //already set - strides will be computed automatically
    explicit
    RangeT(storage_type&& store)
      : store_(std::move(store))
        { 
        computeStrides();
        }

    RangeT&
    operator=(storage_type&& store)
        {
        store_ = std::move(store);
        computeStrides();
        return *this;
        }

    size_type
    extent(size_type i) const { return static_cast<size_type>(store_[i].ext); }

    size_type
    stride(size_type i) const { return store_[i].str; }

    size_type
    r() const { return store_.size(); }

    const value_type&
    operator[](size_type i) const { return store_[i]; }

    value_type&
    operator[](size_type i) { return store_[i]; }

    const value_type&
    at(size_type i) const { return store_.at(i); }

    value_type&
    at(size_type i) { return store_.at(i); }

    size_type
    size() const { return store_.size(); }

    void
    clear() { return store_.clear(); }

    bool
    empty() const { return store_.empty(); }

    extent_type&
    front() { return store_.front().ext; }

    const value_type&
    front() const { return store_.front(); }

    value_type&
    back() { return store_.back(); }

    const value_type&
    back() const { return store_.back(); }

    value_type*
    data() { return store_.data(); }

    const value_type*
    data() const { return store_.data(); }

    const storage_type&
    store() const { return store_; }

    void
    swap(RangeT& other)
        {
        store_.swap(other.store_);
        }

    template<typename Indexable>
    void 
    init(const Indexable& v);

    void 
    init(std::initializer_list<extent_type> il);

    void
    resize(size_type nsize) { store_.resize(nsize); }

    void 
    computeStrides()
        {
        size_type str = 1;
        for(auto& i : store_)
            {
            i.str = str;
            str *= static_cast<size_type>(i.ext);
            }
        }
    };

namespace detail {

template<typename T>
class IListAdapter
    {
    std::initializer_list<T> ilist_;
    public:

    IListAdapter(std::initializer_list<T> il) : ilist_(il) { }

    size_t
    size() const { return ilist_.size(); }

    T&
    operator[](size_t j) { return *(ilist_.begin()+j); }

    const T&
    operator[](size_t j) const { return *(ilist_.begin()+j); }
    };

} //namespace detail

template<typename extent_type>
template<typename Indexable>
void RangeT<extent_type>::
init(const Indexable& v)
    {
    store_.resize(v.size());
    size_type str = 1;
    for(size_type i = 0; i < v.size(); ++i)
        {
        store_[i].ext = v[i];
        store_[i].str = str;
        str *= static_cast<size_type>(v[i]);
        }
    }

template<typename extent_type>
void RangeT<extent_type>::
init(std::initializer_list<extent_type> il)
    {
    init(detail::IListAdapter<extent_type>{il});
    }


template<typename extent_type>
std::ostream&
operator<<(std::ostream& s, const RangeT<extent_type>& r)
    {
    s << "dims: ";
    for(Range::size_type i = 0; i < r.r(); ++i) s << r.dim(i) << " ";
    s << "strs: ";
    for(Range::size_type i = 0; i < r.r(); ++i) s << r.stride(i) << " ";
    return s;
    }

//template<typename extent_type>
//void
//write(std::ostream& s, const ExtStr<extent_type>& xs)
//    {
//    itensor::write(s,xs.ext);
//    itensor::write(s,xs.str);
//    }
//
//template<typename extent_type>
//void
//read(std::istream& s, ExtStr<extent_type>& xs)
//    {
//    itensor::read(s,xs.ext);
//    itensor::read(s,xs.str);
//    }

template<typename extent_type>
void
write(std::ostream& s, const RangeT<extent_type>& r)
    {
    itensor::write(s,r.store());
    }

template<typename extent_type>
void
read(std::istream& s, RangeT<extent_type>& r)
    {
    using storage_type = typename RangeT<extent_type>::storage_type;
    storage_type store;
    itensor::read(s,store);
    r = RangeT<extent_type>(std::move(store));
    }

//class RangeRef
//    {
//    using extent_type = Range::extent_type;
//    using value_type = index;
//    using size_type = Range::size_type;
//    private:
//    value_type* inds_;
//    size_type rank_;
//    public:
//
//    RangeRef(value_type* prange,
//             size_type rank)
//        :
//        inds_(prange),
//        rank_(rank)
//        { }
//
//    size_type
//    dim(size_type i) const { return inds_[i].dim; }
//
//    size_type
//    stride(size_type i) const { return inds_[i].stride; }
//
//    size_type
//    r() const { return rank_; }
//
//    const index&
//    operator[](size_type i) const { return inds_[i]; }
//
//    index&
//    operator[](size_type i) { return inds_[i]; }
//    };

namespace detail {
template<typename RangeT, typename Iterable>
size_t
indIterable(const RangeT& r, const Iterable& inds)
    {
    size_t I  = 0, 
           ri = 0;
    for(auto& ii : inds)
        {
        I += r.stride(ri) * ii;
        ++ri;
        }
    return I;
    }

template<typename RangeT>
struct ComputeInd
    {
    using size_type = typename RangeT::size_type;

    const RangeT& r;
    ComputeInd(const RangeT& r_) : r(r_) { }

    template<typename... Inds>
    size_type
    operator()(size_type first, Inds... rest) const 
        { 
        return ind__<0>(first,rest...);
        }

    private:
    template <size_type i, typename... Inds>
    size_type
    ind__(size_type first, Inds... rest) const
        {
        return first*r.stride(i) + ind__<i+1>(rest...);
        }

    template <size_type i>
    size_type
    ind__(size_type first) const
        {
        return first*r.stride(i);
        }
    };

} //namespace detail

template<typename RangeT, typename U>
size_t
ind(const RangeT& r, const std::vector<U>& inds)
    {
    return detail::indIterable(r,inds);
    }

template<typename RangeT, typename U, size_t size>
size_t
ind(const RangeT& r, const std::array<U,size>& inds)
    {
    return detail::indIterable(r,inds);
    }

template<typename RangeT, typename U, size_t size>
size_t
ind(const RangeT& r, const VarArray<U,size>& inds)
    {
    return detail::indIterable(r,inds);
    }

template<typename RangeT, typename U, size_t size>
size_t
ind(const RangeT& r, const InfArray<U,size>& inds)
    {
    return detail::indIterable(r,inds);
    }

template<typename RangeT, typename U>
size_t
ind(const RangeT& r, const autovector<U>& inds)
    {
    return detail::indIterable(r,inds);
    }

template<typename RangeT>
size_t
ind(const RangeT& r, size_t i0)
    {
    return detail::ComputeInd<RangeT>(r)(i0);
    }

template<typename RangeT, typename... Inds>
size_t
ind(const RangeT& r, size_t i0, size_t i1, Inds... inds)
    {
    return detail::ComputeInd<RangeT>(r)(i0,i1,inds...);
    }

template<typename RangeT>
size_t
area(const RangeT& r)
    { 
    if(r.r()==0) return 1ul;
    //TODO: this won't work if we allow Ranges to be constructed for slicing
    auto last = r.r()-1;
    auto A = r.extent(last)*r.stride(last);
    return A;
    }


} //namespace itensor

#endif
