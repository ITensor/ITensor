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
#include "itensor/tensor/rangeiter.h"

namespace itensor {

template<typename extent_type>
class RangeT;

template<typename extent_type>
class RangeBuilderT;

using Range = RangeT<size_t>;
using RangeBuilder = RangeBuilderT<size_t>;

//Storage type for RangeT
template<typename extent_type>
struct ExtStr
    {
    using size_type = size_t;
    extent_type ext = extent_type{}; //convertible to size_type
    size_type   str = 0; //stride

    ExtStr(extent_type e, size_type s) : ext(e), str(s) { }

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
    using extent_type = extent_type_;
    using value_type = ExtStr<extent_type>;
    using size_type = typename value_type::size_type;
    using storage_type = InfArray<value_type,8ul>;
    using iterator = RangeIter<RangeT>;
    using const_iterator = iterator;
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
    RangeT(std::vector<extent_type> const& v)
        {
        init(v);
        }

    template<size_t size>
    explicit
    RangeT(std::array<extent_type,size> const& a)
        {
        init(a);
        }

    explicit
    RangeT(storage_type && store)
      : store_(std::move(store))
        { }

    RangeT(std::initializer_list<extent_type> ii) { init(ii); }

    //0-indexed
    size_type
    extent(size_type i) const { return static_cast<size_type>(store_[i].ext); }

    //0-indexed
    size_type
    stride(size_type i) const { return store_[i].str; }

    size_type
    r() const { return store_.size(); }

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

    void
    clear() { return store_.clear(); }

    bool
    empty() const { return store_.empty(); }

    value_type*
    data() { return store_.data(); }

    const value_type*
    data() const { return store_.data(); }

    storage_type const&
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

    iterator
    begin() const { return RangeIter<RangeT>(*this); }

    iterator
    end() const { return RangeIter<RangeT>::makeEnd(*this); }

    //Compute strides from extents
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

template<typename extent_type>
template<typename Indexable>
void RangeT<extent_type>::
init(Indexable const& v)
    {
    store_.resize(v.size());
    size_type str = 1;
    for(decltype(v.size()) i = 0; i < v.size(); ++i)
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
    store_.resize(il.size());
    size_type str = 1;
    size_type i = 0;
    for(auto& el : il)
        {
        store_[i].ext = el;
        store_[i].str = str;
        str *= static_cast<size_type>(el);
        ++i;
        }
    }


template<typename extent_type>
std::ostream&
operator<<(std::ostream& s, RangeT<extent_type> const& r)
    {
    s << "exts: ";
    for(decltype(r.r()) i = 0; i < r.r(); ++i) s << r.extent(i) << " ";
    s << "strs: ";
    for(decltype(r.r()) i = 0; i < r.r(); ++i) s << r.stride(i) << " ";
    return s;
    }


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


template<typename extent_type_>
class RangeBuilderT
    {
    public:
    using extent_type = extent_type_;
    using range_type = RangeT<extent_type>;
    using size_type = typename range_type::size_type;
    using storage_type = typename range_type::storage_type;
    using value_type = typename storage_type::value_type;
    private:
    storage_type store_;
    bool auto_compute_strides_ = true;
    public:

    RangeBuilderT() { }

    RangeBuilderT(size_type size)
      : store_(size)
        { }

    explicit
    operator range_type()
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
    setExtent(size_type j, extent_type const& e)
        {
        store_[j].ext = e;
        }

    void
    setStride(size_type j, size_type s)
        {
        auto_compute_strides_ = false;
        store_[j].str = s;
        }
    };



namespace detail {

    template<typename RangeT, typename Iterable>
    auto
    indIterable(RangeT const& r, Iterable const& inds)
        {
        using size_type = typename RangeT::size_type;
        size_type I  = 0, 
                  ri = 0;
        for(auto& ii : inds)
            {
#ifdef DEBUG
            if(ri >= size_type(r.r()))
                Error("Container-Range mismatch in ind(...)");
#endif
            I += r.stride(ri) * ii;
            ++ri;
            }
        return I;
        }

    template<typename RangeT>
    struct ComputeInd
        {
        using size_type = typename RangeT::size_type;

        RangeT const& r;
        ComputeInd(RangeT const& r_) : r(r_) { }

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

template<typename RangeT, typename Iterable>
auto
ind(RangeT const& r, Iterable const& inds)
    {
    return detail::indIterable(r,inds);
    }

template<typename RangeT, typename U>
auto
ind(RangeT const& r, std::vector<U> const& inds)
    {
    return detail::indIterable(r,inds);
    }

template<typename RangeT, typename U, size_t size>
auto
ind(RangeT const& r, std::array<U,size> const& inds)
    {
    return detail::indIterable(r,inds);
    }

template<typename RangeT, typename U, size_t size>
auto
ind(RangeT const& r, VarArray<U,size> const& inds)
    {
    return detail::indIterable(r,inds);
    }

template<typename RangeT, typename U, size_t size>
auto
ind(RangeT const& r, InfArray<U,size> const& inds)
    {
    return detail::indIterable(r,inds);
    }

template<typename RangeT, typename U>
auto
ind(RangeT const& r, autovector<U> const& inds)
    {
    return detail::indIterable(r,inds);
    }

template<typename RangeT>
auto
ind(RangeT const& r, size_t i0)
    {
    return detail::ComputeInd<RangeT>(r)(i0);
    }

template<typename RangeT, typename... Inds>
auto
ind(RangeT const& r, size_t i0, size_t i1, Inds... inds)
    {
    return detail::ComputeInd<RangeT>(r)(i0,i1,inds...);
    }

template<typename RangeT>
auto
area(RangeT const& R)
    { 
    using size_type = typename RangeT::size_type;
    size_type A = 1;
    for(decltype(R.r()) n = 0; n < R.r(); ++n)
        {
        A *= R.extent(n);
        }
    return A;
    }

//
//A range R is contiguous if collecting
//all possible outputs of ind(R,...) yields
//the set {0,1,...,area(R)-1} (though 
//in no particular order)
//For this to be true, sufficient that
//the max possible output of ind(R,...)
//equals (area(R)-1)
//
//Proof:
// Ranges always start at ind(R,{0,0,0...})=0;
// area(R) gives the number of outputs;
// IF the max output is area(R)-1 and
// there are area(R) outputs, the only
// set fulfilling this is {0,1,...,area(R)-1})
//
template<typename RangeT>
bool
isContiguous(RangeT const& R)
    {
    using size_type = typename RangeT::size_type;
    size_type max_ind = 0,
              area = 1;
    for(decltype(R.r()) n = 0; n < R.r(); ++n)
        {
        max_ind += R.stride(n)*(R.extent(n)-1);
        area *= R.extent(n);
        }
    return (1+max_ind) == area;
    }


} //namespace itensor

#endif
