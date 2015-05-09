//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_RANGE_H
#define __ITENSOR_RANGE_H

#include "autovector.h"
#include <iostream>
#include <array>

namespace itensor {

class Range
    {
    public:
    struct index 
        { 
        long dim = 0, 
             stride = 0; 
        index() { }
        index(long dim_, long stride_) : dim(dim_),stride(stride_) { }
        explicit operator long() const { return dim; }
        };
    using value_type = index;
    using storage_type = std::vector<index>;

    Range() { }

    template <typename... Dims>
    Range(long d0, Dims... rest)
        {
        std::array<long,1+sizeof...(rest)> dims = {{d0,static_cast<long>(rest)...}};
        init(dims);
        }

    template <typename I>
    Range(const std::vector<I>& v)
        {
        init(v);
        }

    template <typename I, size_t size>
    Range(const std::array<I,size>& a)
        {
        init(a);
        }

    long
    dim(long i) const { return inds_[i].dim; }
    long
    stride(long i) const { return inds_[i].stride; }
    long
    r() const { return inds_.size(); }

    const index&
    operator[](long i) const { return inds_[i]; }
    index&
    operator[](long i) { return inds_[i]; }

    const index&
    at(long i) const { return inds_.at(i); }
    index&
    at(long i) { return inds_.at(i); }

    size_t
    size() const { return inds_.size(); }

    void
    clear() { return inds_.clear(); }

    bool
    empty() const { return inds_.empty(); }

    index&
    front() { return inds_.front(); }
    const index&
    front() const { return inds_.front(); }

    index&
    back() { return inds_.back(); }
    const index&
    back() const { return inds_.back(); }

    void
    swap(Range& other)
        {
        inds_.swap(other.inds_);
        }

    template<typename Indexable>
    void 
    init(const Indexable& v)
        {
        inds_.resize(v.size());
        long len = 1;
        for(size_t i = 0; i < size_t(v.size()); ++i)
            {
            inds_[i] = index(long{v[i]},len);
            len *= long{v[i]};
            }
        }

    private:

    storage_type inds_;
    };

inline
std::ostream&
operator<<(std::ostream& s, const Range& r)
    {
    s << "dim: ";
    for(long i = 0; i < r.r(); ++i) s << r.dim(i) << " ";
    s << "stride: ";
    for(long i = 0; i < r.r(); ++i) s << r.stride(i) << " ";
    return s;
    }

class RangeRef
    {
    public:
    using index = Range::index;
    using value_type = index;

    RangeRef(value_type* prange,
             long rank)
        :
        inds_(prange),
        rank_(rank)
        { }

    long
    dim(long i) const { return inds_[i].dim; }
    long
    stride(long i) const { return inds_[i].stride; }
    long
    r() const { return rank_; }

    const index&
    operator[](long i) const { return inds_[i]; }
    index&
    operator[](long i) { return inds_[i]; }

    private:
    /////////
    value_type* inds_;
    long rank_;
    /////////
    };

namespace detail {
template<typename RangeT, typename Iterable>
long
indIterable(const RangeT& r, const Iterable& inds)
    {
    long ii = 0, 
         i = 0;
    for(auto& j : inds)
        {
        ii += r.stride(i) * j;
        ++i;
        }
    return ii;
    }

template<typename RangeT>
struct ComputeInd
    {
    const RangeT& r;
    ComputeInd(const RangeT& r_) : r(r_) { }

    template<typename... Inds>
    long
    operator()(long first, Inds... rest) const 
        { 
        return ind__<0>(first,rest...);
        }

    private:

    template <long i, typename... Inds>
    long
    ind__(long first, Inds... rest) const
        {
        return first*r.stride(i) + ind__<i+1>(rest...);
        }

    template <long i>
    long
    ind__(long first) const
        {
        return first*r.stride(i);
        }
    };

} //namespace detail

template<typename RangeT, typename U>
long
ind(const RangeT& r, const std::vector<U>& inds)
    {
    return detail::indIterable(r,inds);
    }

template<typename RangeT, typename U, size_t size>
long
ind(const RangeT& r, const std::array<U,size>& inds)
    {
    return detail::indIterable(r,inds);
    }

template<typename RangeT, typename U>
long
ind(const RangeT& r, const autovector<U>& inds)
    {
    return detail::indIterable(r,inds);
    }

template<typename RangeT>
long
ind(const RangeT& r, long i0)
    {
    return detail::ComputeInd<RangeT>(r)(i0);
    }

template<typename RangeT, typename... Inds>
long
ind(const RangeT& r, long i0, long i1, Inds... inds)
    {
    return detail::ComputeInd<RangeT>(r)(i0,i1,inds...);
    }

template<typename RangeT>
long
area(const RangeT& r)
    { 
    if(r.empty()) return 1l;
    //TODO: this won't work if we allow Ranges to be constructed for slicing
    auto last = r.r()-1;
    auto A = r.dim(last)*r.stride(last);
    return A;
    }


} //namespace itensor

#endif
