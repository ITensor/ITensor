//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MATRIX_MATRANGE_H_
#define __ITENSOR_MATRIX_MATRANGE_H_

#include "itensor/tensor/vec.h"

namespace itensor {

class MatRangeIter;

template<size_t start>
struct MatRangeT;

using MatRange = MatRangeT<0ul>;
using MatRange1 = MatRangeT<1ul>;

struct MatRangeType : public RangeType { };

template<size_t start_>
struct MatRangeT : public MatRangeType
    {
    public:
    using size_type = size_t;
    using iterator = MatRangeIter;
    using const_iterator = MatRangeIter;
    size_type rn = 0, //number of rows
              rs = 0, //row stride
              cn = 0, //number of cols
              cs = 0; //column stride

    MatRangeT() { } 

    MatRangeT(size_type rn_, 
             size_type cn_)
        : rn(rn_),rs(1),cn(cn_),cs(rn_) 
        { }

    MatRangeT(size_type rn_, 
             size_type rs_,
             size_type cn_, 
             size_type cs_) 
        : rn(rn_),rs(rs_),cn(cn_),cs(cs_) 
        { }

    size_type constexpr
    start() const { return start_; }

    size_type
    extent(size_type i) const
        {
        return (i==0 ? rn : cn);
        }

    size_type
    stride(size_type i) const
        {
        return (i==0 ? rs : cs);
        }

    size_type
    r() const { return 2; }

    iterator
    begin() const;

    iterator
    end() const;

    void
    clear()
        {
        rn = 0; rs = 0; cn = 0; cs = 0;
        }

    };

template<size_t S>
size_t
rank(MatRangeT<S> const& R) { return 2ul; }

//make MatRange with same extents 
//but usual strides
template<size_t S>
MatRangeT<S>
normalRange(MatRangeT<S> const& mr)
    {
    return MatRangeT<S>{mr.rn,mr.cn};
    }

template<size_t S>
MatRangeT<S>
transpose(MatRangeT<S> const& mr)
    {
    return MatRangeT<S>{mr.cn,mr.cs,mr.rn,mr.rs};
    }

//template<size_t S>
//auto
//offset(MatRangeT<S> const& mr, 
//       size_t i1,
//       size_t i2)
//    -> typename MatRangeT<S>::size_type
//    {
//    return i1*mr.rs+i2*mr.cs;
//    }
//
//template<size_t S, typename Iterable>
//auto
//offset(MatRangeT<S> const& mr, Iterable const& inds)
//    -> stdx::if_compiles_return<typename MatRangeT<S>::size_type,decltype(inds[0])>
//    {
//    assert(inds.size()==2);
//    return offset(mr,inds[0],inds[1]);
//    }

template<size_t S>
auto
area(MatRangeT<S> const& mr)
    -> typename MatRangeT<S>::size_type
    {
    return mr.rn * mr.cn;
    }

template<size_t S>
bool
operator==(MatRangeT<S> const& a, MatRangeT<S> const& b)
    {
    return (a.rn==b.rn && a.rs==b.rs && a.cn==b.cn && a.cs==b.cs);
    }

template<size_t S>
bool
operator!=(MatRangeT<S> const& a, MatRangeT<S> const& b) { return !operator==(a,b); }

template<size_t S>
bool
isTransposed(MatRangeT<S> const& i) { return (i.rs==i.cn && i.cs==1); }

template<size_t S>
bool
isNormal(MatRangeT<S> const& i) { return (i.rs==1 && i.cs==i.rn); }

template<size_t S>
bool
isContiguous(MatRangeT<S> const& i) { return isNormal(i) || isTransposed(i); }

template<size_t S>
std::ostream&
operator<<(std::ostream& s, MatRangeT<S> const& mr)
    {
    s << "(rn="<< mr.rn <<",rs="<< mr.rs <<",cn="<< mr.cn <<",cs="<< mr.cs << ")";
    return s;
    }


class MatRangeIter
    { 
    public:
    using size_type = size_t;
    using offset_type = size_type;
    using ind_type = offset_type;
    using iterator_category = std::forward_iterator_tag;
    using range_type = MatRange;
    private:
    size_type count_ = 0;
    offset_type off_ = 0;
    range_type range_; 
    public: 

    MatRangeIter() { }

    explicit
    MatRangeIter(range_type const& mr) : range_(mr) { }  

    offset_type
    offset() const { return off_; }

    MatRangeIter& 
    operator++() { increment(); return *this; } 

    MatRangeIter 
    operator++(int) { auto prev = *this; increment(); return prev; } 

    MatRangeIter const& 
    operator*() { return *this; }  

    bool
    operator==(MatRangeIter const& o) const { return count_ == o.count_; }

    bool
    operator!=(MatRangeIter const& o) const { return !operator==(o); }

    private:

    void
    increment()
        {
        //print("off_ ",off_," -> ");
        ++count_;
        off_ += range_.rs;
        //println(off_);
        if((count_%range_.rn) == 0)
            {
            //printfln("count_ mod %d == 0",range_.rn);
            off_ -= range_.rs*range_.rn;
            off_ += range_.cs;
            //println("now off_ = ",off_);
            }
        }

    public:
    //For developer use only; for making end iterator
    MatRangeIter static
    makeEnd(range_type const& r)
        {
        MatRangeIter end;
        end.count_ = r.rn*r.cn;
        //end.off_ = 1 + r.rs*(r.rn-1) + r.cs*(r.cn-1);
        //println("end.off_ = ",end.off_);
        return end;
        }
    }; 

template<size_t S>
auto MatRangeT<S>::
begin() const -> iterator { return iterator(*this); }

template<size_t S>
auto MatRangeT<S>::
end() const -> iterator { return iterator::makeEnd(*this); }

} //namespace itensor

#endif
