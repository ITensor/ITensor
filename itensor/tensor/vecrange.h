//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MATRIX_VECRANGE_H_
#define __ITENSOR_MATRIX_VECRANGE_H_

#include <iterator> 
#include "itensor/util/print.h"
#include "itensor/util/error.h"
#include "itensor/tensor/range.h"

namespace itensor {

class VecRangeIter;


template<size_t start>
class VecRangeT;

using VecRange = VecRangeT<0ul>;
using VecRange1 = VecRangeT<1ul>;

struct VecRangeType : public RangeType { };

template<size_t start_>
class VecRangeT : public VecRangeType
    {
    public:
    using size_type = size_t;
    using iterator = VecRangeIter;
    using const_iterator = VecRangeIter;
    private:
    size_type ext_ = 0,
              stride_ = 1;
    public:

    VecRangeT() { } 

    explicit
    VecRangeT(size_type extent)
      : ext_(extent),
        stride_(1)
        { }

    VecRangeT(size_type extent,
             size_type stride)
      : ext_(extent),
        stride_(stride)
        { }

    size_type constexpr
    start() const { return start_; }

    size_type
    extent() const { return ext_; }

    size_type
    stride() const { return stride_; }

    size_type
    extent(size_type i) const 
        { 
#ifdef DEBUG
        if(i != 0) Error("i out of range in VecRangeT::extent(i)");
#endif
        return ext_; 
        }

    size_type
    stride(size_type i) const 
        {
#ifdef DEBUG
        if(i != 0) Error("i out of range in VecRangeT::stride(i)");
#endif
        return stride_; 
        }

    size_type
    r() const { return 1; }

    size_type
    size() const { return 1; }

    iterator
    begin() const;

    iterator
    end() const;

    void
    clear()
        {
        ext_ = 0;
        stride_ = 1;
        }

    void
    read(std::istream& s)
        {
        itensor::read(s,ext_);
        itensor::read(s,stride_);
        }

    void
    write(std::ostream& s) const
        {
        itensor::write(s,ext_);
        itensor::write(s,stride_);
        }

    };

template<size_t S>
size_t
rank(VecRangeT<S> const& R) { return 1ul; }

//template<size_t S>
//size_t
//order(VecRangeT<S> const& R) { return 1ul; }

//make VecRange with same extent but stride()==1
template<size_t S>
VecRangeT<S>
normalRange(VecRangeT<S> const& vr)
    {
    return VecRangeT<S>{vr.extent()};
    }

////0-indexed
//auto inline
//offset(VecRange const& vr, VecRange::size_type ind)
//    -> VecRange::size_type
//    {
//    return vr.stride()*ind;
//    }
//
//template<typename Inds>
//auto
//offset(VecRange const& vr, Inds const& i)
//    -> stdx::if_compiles_return<VecRange::size_type,decltype(i[0])>
//    {
//    return vr.stride()*i[0];
//    }

template<size_t S>
auto
area(VecRangeT<S> const& vr)
    -> decltype(vr.extent())
    {
    return vr.extent();
    }

template<size_t S>
bool
isNormal(VecRangeT<S> const& vr) { return vr.stride()==1; }

template<size_t S>
bool
isContiguous(VecRangeT<S> const& vr) { return vr.stride()==1; }

template<size_t S>
std::ostream&
operator<<(std::ostream& s, VecRangeT<S> const& vr)
    {
    s << "(extent="<< vr.extent() <<",stride="<< vr.stride() << ")";
    return s;
    }


class VecRangeIter
    { 
    public:
    using size_type = size_t;
    using offset_type = size_type;
    using ind_type = size_type;
    using iterator_category = std::forward_iterator_tag;
    using range_type = VecRange;
    private:
    offset_type off_ = 0;
    range_type range_;
    public: 

    VecRangeIter() { }

    explicit
    VecRangeIter(range_type const& vr) : range_(vr) { }  

    offset_type
    offset() const { return off_; }

    ind_type
    index() const { return range_.start()+off_; }

    size_type
    stride() const { return range_.stride(); }

    VecRangeIter& 
    operator++() { off_ += stride(); return *this; } 

    VecRangeIter 
    operator++(int) { auto old = *this; off_ += stride(); return old; } 

    VecRangeIter& 
    operator+=(size_type x) { off_ += x * stride(); return *this; } 

    VecRangeIter& 
    operator--() { off_ -= stride(); return *this; } 

    VecRangeIter 
    operator--(int) { auto old = *this; off_ -= stride(); return old; } 

    VecRangeIter& 
    operator-=(size_type x) { off_ -= x * stride(); return *this; } 

    VecRangeIter const& 
    operator*() { return *this; }  

    VecRangeIter static
    makeEnd(range_type const& r)
        {
        auto end = VecRangeIter(r);
        end.off_ = r.stride()*r.extent();
        return end;
        }
    }; 

bool inline
operator==(VecRangeIter const& x, VecRangeIter const& y) 
    { assert(x.stride() == y.stride()); return x.offset() == y.offset(); } 
bool inline
operator!=(VecRangeIter const& x, VecRangeIter const& y) 
    { assert(x.stride() == y.stride()); return x.offset() != y.offset(); } 
bool inline
operator<(VecRangeIter const& x, VecRangeIter const& y) 
    { assert(x.stride() == y.stride()); return x.offset() < y.offset(); } 

template<size_t S>
auto VecRangeT<S>::
begin() const -> iterator { return iterator(*this); }

template<size_t S>
auto VecRangeT<S>::
end() const -> iterator { return iterator::makeEnd(*this); }

} //namespace itensor

#endif
