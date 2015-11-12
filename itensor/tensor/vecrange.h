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

class VecRange : public RangeType
    {
    public:
    using size_type = size_t;
    using iterator = VecRangeIter;
    using const_iterator = VecRangeIter;
    private:
    size_type ext_ = 0,
              stride_ = 1;
    public:

    VecRange() { } 

    explicit
    VecRange(size_type extent)
      : ext_(extent),
        stride_(1)
        { }

    VecRange(size_type extent,
             size_type stride)
      : ext_(extent),
        stride_(stride)
        { }

    size_type
    extent() const { return ext_; }

    size_type
    stride() const { return stride_; }

    size_type
    extent(size_type i) const 
        { 
#ifdef DEBUG
        if(i != 0) Error("i out of range in VecRange::extent(i)");
#endif
        return ext_; 
        }

    size_type
    stride(size_type i) const 
        {
#ifdef DEBUG
        if(i != 0) Error("i out of range in VecRange::stride(i)");
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

    };

VecRange::size_type inline
rank(VecRange const& R) { return 1; }

//make VecRange with same extent but stride()==1
VecRange inline
normalRange(VecRange const& vr)
    {
    return VecRange{vr.extent()};
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

auto inline
area(VecRange const& vr)
    -> VecRange::size_type
    {
    return vr.extent();
    }

bool inline
isNormal(VecRange const& vr) { return vr.stride()==1; }

bool inline
isContiguous(VecRange const& vr) { return vr.stride()==1; }

inline std::ostream&
operator<<(std::ostream& s, VecRange const& vr)
    {
    s << "(extent="<< vr.extent() <<",stride="<< vr.stride() << ")";
    return s;
    }


class VecRangeIter
    { 
    public:
    using size_type = size_t;
    using offset_type = size_type;
    using iterator_category = std::forward_iterator_tag;
    using range_type = VecRange;
    private:
    offset_type off_ = 0;
    size_type stride_ = 0; 
    public: 

    VecRangeIter() { }

    explicit
    VecRangeIter(size_type stride) : stride_(stride) { }  

    explicit
    VecRangeIter(range_type const& vr) : stride_(vr.stride()) { }  

    offset_type
    offset() const { return off_; }

    size_type
    stride() const { return stride_; }

    VecRangeIter& 
    operator++() { off_ += stride_; return *this; } 

    VecRangeIter 
    operator++(int) { auto old = *this; off_ += stride_; return old; } 

    VecRangeIter& 
    operator+=(size_type x) { off_ += x * stride_; return *this; } 

    VecRangeIter& 
    operator--() { off_ -= stride_; return *this; } 

    VecRangeIter 
    operator--(int) { auto old = *this; off_ -= stride_; return old; } 

    VecRangeIter& 
    operator-=(size_type x) { off_ -= x * stride_; return *this; } 

    VecRangeIter const& 
    operator*() { return *this; }  

    VecRangeIter static
    makeEnd(range_type const& r)
        {
        VecRangeIter end;
        end.off_ = r.stride()*r.extent();
        end.stride_ = r.stride();
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

VecRange::iterator inline VecRange::
begin() const { return iterator(*this); }

VecRange::iterator inline VecRange::
end() const { return iterator::makeEnd(*this); }

} //namespace itensor

#endif
