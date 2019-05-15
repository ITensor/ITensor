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
#ifndef __ITENSOR_RANGE_H
#define __ITENSOR_RANGE_H

#include "itensor/util/stdx.h"
#include "itensor/util/readwrite.h"
#include "itensor/tensor/rangeiter.h"

namespace itensor {

//Parent type for all ranges
class RangeType { };

template<typename index_type,
         size_t start = 0>
class RangeT;

template<typename range_type>
class RangeBuilderT;

//0-indexed
using Range = RangeT<size_t,0ul>;
using RangeBuilder = RangeBuilderT<Range>;

//1-indexed
using Range1 = RangeT<size_t,1ul>;
using Range1Builder = RangeBuilderT<Range1>;

template<typename Derived>
struct isRange
    {
    bool static constexpr value = std::is_base_of<RangeType,Derived>::value;
    constexpr operator bool() const noexcept { return value; }
    };

//Storage type for RangeT
template<typename index_type>
struct IndStr
    {
    using size_type = size_t;
    index_type  ind = index_type{}; //convertible to size_type
    size_type   str = 0; //stride

    IndStr(index_type i, size_type s) : ind(i), str(s) { }

    IndStr() { }

    void
    write(std::ostream& s) const
        {
        itensor::write(s,ind);
        itensor::write(s,str);
        }
    void
    read(std::istream& s)
        {
        itensor::read(s,ind);
        itensor::read(s,str);
        }
    };

template<typename index_type_, size_t start_>
class RangeT : public RangeType
    {
    public:
    using index_type = index_type_;
    using value_type = IndStr<index_type>;
    using size_type = typename value_type::size_type;
    using storage_type = InfArray<value_type,11ul>;
    using iterator = RangeIter<RangeT>;
    using const_iterator = iterator;
    private:
    storage_type store_;
    public:

    RangeT() { }

    template <typename... Inds>
    explicit
    RangeT(index_type i0, 
           Inds&&... rest)
        {
        std::array<index_type,1+sizeof...(rest)> inds = 
            {{i0,static_cast<index_type>(rest)...}};
        init(inds);
        }

    explicit
    RangeT(storage_type && store)
      : store_(std::move(store))
        { }

    explicit
    RangeT(std::vector<index_type> const& inds)
        {
        init(inds);
        }

    template<size_t N>
    explicit
    RangeT(std::array<index_type,N> const& inds)
        {
        init(inds);
        }

    RangeT(std::initializer_list<index_type> ii) { init(ii); }

    size_t constexpr
    start() const { return start_; }

    size_type
    extent(size_type i) const { return static_cast<size_type>(store_[i].ind); }

    size_type
    stride(size_type i) const { return store_[i].str; }

    index_type const&
    index(size_type i) const { return store_[i].ind; }

    index_type &
    index(size_type i) { return store_[i].ind; }

    size_type
    order() const { return store_.size(); }

    // Deprecated
    size_type
    r() const { return this->order(); }

    value_type const&
    operator[](size_type i) const { return store_[i]; }

    value_type &
    operator[](size_type i) { return store_[i]; }

    value_type const&
    at(size_type i) const { return store_.at(i); }

    value_type &
    at(size_type i) { return store_.at(i); }

    size_type
    size() const { return store_.size(); }

    value_type const&
    front() const { return store_.front(); }

    value_type const&
    back() const { return store_.back(); }

    void
    clear() { return store_.clear(); }

    bool
    empty() const { return store_.empty(); }

    value_type*
    data() { return store_.data(); }

    value_type const*
    data() const { return store_.data(); }

    storage_type const&
    store() const { return store_; }

    void
    swap(RangeT& other) { store_.swap(other.store_); }

    template<typename Container>
    void 
    init(Container const& c);

    void
    resize(size_type nsize) { store_.resize(nsize); }

    iterator
    begin() const { return iterator(*this); }

    iterator
    end() const { return iterator::makeEnd(*this); }

    //Compute strides from extents
    void 
    computeStrides()
        {
        size_type str = 1;
        for(auto& i : store_)
            {
            i.str = str;
            str *= static_cast<size_type>(i.ind);
            }
        }
    };

namespace detail {

template<typename range_type,
         typename Indexable,
         typename storage_type = typename range_type::storage_type>
auto
initImpl(stdx::choice<2>,
         Indexable && v,
         storage_type & store_)
    -> stdx::if_compiles_return<void,decltype(v[0]),decltype(v.size())>
    {
    using size_type = typename range_type::size_type;
    store_.resize(v.size());
    size_type str = 1;
    for(decltype(v.size()) i = 0; i < v.size(); ++i)
        {
        store_[i].ind = v[i];
        store_[i].str = str;
        str *= static_cast<size_type>(v[i]);
        }
    }

template<typename range_type,
         typename Iterable,
         typename storage_type = typename range_type::storage_type>
auto
initImpl(stdx::choice<1>,
         Iterable && v,
         storage_type & store_)
    -> stdx::if_compiles_return<void,decltype(v.begin())>
    {
    using size_type = typename range_type::size_type;
    store_.resize(v.size());
    size_type str = 1;
    size_type i = 0;
    for(auto& vel : v)
        {
        store_[i].ind = vel;
        store_[i].str = str;
        str *= static_cast<size_type>(vel);
        ++i;
        }
    }

} //namespace detail

template<typename index_type,size_t start>
template<typename Container>
void RangeT<index_type,start>::
init(Container const& c)
    {
    detail::initImpl<RangeT>(stdx::select_overload{},c,store_);
    }

template<typename index_type, size_t start>
auto
order(RangeT<index_type,start> const& R) -> decltype(R.size()) { return R.size(); }

//template<typename index_type, size_t start>
//auto
//order(RangeT<index_type,start> const& R) -> decltype(R.size()) { return R.size(); }


template<typename index_type,size_t start>
std::ostream&
operator<<(std::ostream& s, RangeT<index_type,start> const& r)
    {
    s << "exts: ";
    for(decltype(r.order()) i = 0; i < r.order(); ++i) s << r.extent(i) << " ";
    s << "strs: ";
    for(decltype(r.order()) i = 0; i < r.order(); ++i) s << r.stride(i) << " ";
    return s;
    }


template<typename I, size_t S>
void
write(std::ostream& s, RangeT<I,S> const& r)
    {
    itensor::write(s,r.store());
    }

template<typename I, size_t S>
void
read(std::istream& s, RangeT<I,S>& r)
    {
    using storage_type = typename RangeT<I,S>::storage_type;
    storage_type store;
    itensor::read(s,store);
    r = RangeT<I,S>(std::move(store));
    }


template<typename range_type_>
class RangeBuilderT
    {
    public:
    using range_type = typename std::decay<range_type_>::type;
    using index_type = typename range_type::index_type;
    using size_type = typename range_type::size_type;
    using storage_type = typename range_type::storage_type;
    using value_type = typename storage_type::value_type;
    private:
    storage_type store_;
    size_t n_ = 0;
    bool auto_compute_strides_ = true;
    public:

    RangeBuilderT(size_type size)
      : store_(size)
        { }

    //"build" and return the range
    //store_ will be moved into the newly created range
    range_type
    build()
        {
        auto res = range_type{std::move(store_)};
        if(auto_compute_strides_) res.computeStrides();
        return res;
        }

    operator bool() const { return !store_.empty(); }

    size_type
    size() const { return store_.size(); }

    void
    resize(size_type newsize) { store_.resize(newsize); }

    size_type const&
    extent(size_type j) const { return static_cast<size_type>(store_[j].ind); }

    size_type
    stride(size_type j) const { return store_[j].str; }

    storage_type &
    store() { return store_; }

    storage_type const&
    store() const { return store_; }

    void
    nextIndex(index_type const& ind)
        {
        store_[n_].ind = ind;
        ++n_;
        }

    void
    nextIndex(index_type && ind)
        {
        store_[n_].ind = std::move(ind);
        ++n_;
        }

    void
    nextIndStr(index_type const& ind,
               size_type str)
        {
        store_[n_].ind = ind;
        store_[n_].str = str;
        auto_compute_strides_ = false;
        ++n_;
        }

    void
    setIndex(size_type j, index_type const& i)
        {
        store_[j].ind = i;
        }

    void
    setStride(size_type j, size_type s)
        {
        store_[j].str = s;
        auto_compute_strides_ = false;
        }

    void
    setIndStr(size_type j, 
              index_type const& ind,
              size_type str)
        {
        store_[j].ind = ind;
        store_[j].str = str;
        auto_compute_strides_ = false;
        }

    void
    sortByIndex()
        {
        auto comp = [](value_type const& i1, value_type const& i2) { return i1.ind > i2.ind; };
        std::sort(store_.begin(),store_.end(),comp);
        auto_compute_strides_ = true;
        }
    };

namespace detail {

    //Implementation of offset for Iterables supporting iteration (begin/end)
    template<typename Range_, typename Iterable>
    auto
    offsetImpl(stdx::choice<1>, Range_ const& r, Iterable const& inds)
        //Constrain this template to only work for inds that have a begin() method
        -> stdx::if_compiles_return<decltype(r.extent(0)),decltype(inds.begin())>
        //...if so make the return type to be decltype(r.extent(1))
        {
        using size_type = typename Range_::size_type;
        auto start = r.start();
        size_type I  = 0, 
                  ri = 0;
        for(auto& ii : inds)
            {
#ifdef DEBUG
            if(ri >= size_type(r.order()))
                Error("Container-Range size mismatch in offset(...)");
#endif
            I += r.stride(ri)*(ii-start);
            ++ri;
            }
        return I;
        }

    //Implementation of offset for Iterables supporting .operator[] and .size()
    template<typename Range_, typename Iterable>
    auto
    offsetImpl(stdx::choice<2>, Range_ const& r, Iterable const& inds)
        //Constrain this template to only work for inds that have .operator[] and .size()
        -> stdx::if_compiles_return<decltype(r.extent(0)),decltype(inds[0]),decltype(inds.size())>
        //...if so make the return type to be decltype(r.extent(0))
        {
        using size_type = decltype(r.extent(0));
        auto start = r.start();
        size_type I  = 0;
        for(decltype(inds.size()) n = 0; n < inds.size(); ++n)
            {
#ifdef DEBUG
            if(static_cast<size_type>(n) >= r.order())
                Error("Container-Range size mismatch in offset(...)");
#endif
            I += r.stride(n)*(inds[n]-start);
            }
        return I;
        }

    //Implementation of offset for Iterables exposing a .offset() method
    template<typename Range_, typename Iterable>
    auto
    offsetImpl(stdx::choice<3>, Range_ const& r, Iterable const& I)
        -> decltype(I.offset())
        {
        return (I.offset()-r.start());
        }

} //namespace detail


template<typename Range_, typename Iterable,
         class=stdx::enable_if_t<isRange<Range_>::value 
                              && not std::is_integral<Iterable>::value>>
auto
offset(Range_ const& r, Iterable const& inds)
    -> typename Range_::size_type
    {
    return detail::offsetImpl(stdx::select_overload{},r,inds);
    }

template<typename Range_, 
         class=stdx::require<isRange<Range_>>,
         typename... Inds>
auto
offset(Range_ const& r, size_t i1, Inds... inds)
    -> decltype(r.stride(0))
    {
#ifdef DEBUG
    if(1+sizeof...(inds) != order(r)) 
        throw std::runtime_error(format("Wrong number of indices passed to TenRef (expected %d got %d)",order(r),1+sizeof...(inds)));
#endif
    auto ia = stdx::make_array(i1,inds...);
    return detail::offsetImpl(stdx::select_overload{},r,ia);
    }

template<typename I, size_t S>
auto
dim(RangeT<I,S> const& R)
    -> decltype(R.extent(0))
    { 
    using size_type = decltype(R.size());
    size_type A = 1;
    for(decltype(R.order()) n = 0; n < R.order(); ++n)
        {
        A *= R.extent(n);
        }
    return A;
    }

//make Range with same extents but
//normal (unsliced) strides
template<typename I, size_t S>
Range
normalRange(RangeT<I,S> const& R)
    {
    auto rb = RangeBuilder(R.order());
    for(decltype(R.order()) n = 0; n < R.order(); ++n)
        rb.nextIndex(R.extent(n));
    return rb.build();
    }


//
//A range R is contiguous if collecting
//all possible outputs of offset(R,...) yields
//the set {0,1,...,dim(R)-1} (though 
//in no particular order)
//For this to be true, sufficient that
//the max possible output of offset(R,...)
//equals (dim(R)-1)
//
//Proof:
// Ranges always start at offset(R,{1,1,1...})=0;
// dim(R) gives the number of outputs;
// IF the max output is dim(R)-1 and
// there are dim(R) outputs, the only
// set fulfilling this is {0,1,...,dim(R)-1})
//
template<typename I, size_t S>
bool
isContiguous(RangeT<I,S> const& R)
    {
    using size_type = decltype(R.size());
    size_type max_offset = 0,
              dim = 1;
    for(decltype(R.order()) n = 0; n < R.order(); ++n)
        {
        max_offset += R.stride(n)*(R.extent(n)-1);
        dim *= R.extent(n);
        }
    return (1+max_offset) == dim;
    }


} //namespace itensor

#endif
