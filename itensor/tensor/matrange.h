//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MATRIX_MATRANGE_H_
#define __ITENSOR_MATRIX_MATRANGE_H_

#include "itensor/tensor/vec.h"

namespace itensor {

class MatRangeIter;

struct MatRange
    {
    public:
    using size_type = size_t;
    using iterator = MatRangeIter;
    using const_iterator = MatRangeIter;
    size_type rn = 0, //number of rows
              rs = 0, //row stride
              cn = 0, //number of cols
              cs = 0; //column stride

    MatRange() { } 

    MatRange(size_type rn_, 
             size_type cn_)
        : rn(rn_),rs(1),cn(cn_),cs(rn_) 
        { }

    MatRange(size_type rn_, 
             size_type rs_,
             size_type cn_, 
             size_type cs_) 
        : rn(rn_),rs(rs_),cn(cn_),cs(cs_) 
        { }

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

    };

//make MatRange with same extents 
//but usual strides
MatRange inline
normalRange(MatRange const& mr)
    {
    return MatRange{mr.rn,mr.cn};
    }

MatRange inline
transpose(MatRange const& mr)
    {
    return MatRange{mr.cn,mr.cs,mr.rn,mr.rs};
    }

//0-indexed
auto inline
offset(MatRange const& mr, 
       size_t i1,
       size_t i2)
    {
    return i1*mr.rs+i2*mr.cs;
    }

//0-indexed
template<typename Iterable>
auto
offset(MatRange const& mr, Iterable const& inds)
    -> stdx::if_compiles_return<decltype(inds[0]),MatRange::size_type>
    {
    assert(inds.size()==2);
    return offset(mr,inds[0],inds[1]);
    }

auto inline
area(MatRange const& mr)
    {
    return mr.rn * mr.cn;
    }

bool inline
operator==(MatRange const& a, MatRange const& b)
    {
    return (a.rn==b.rn && a.rs==b.rs && a.cn==b.cn && a.cs==b.cs);
    }

bool inline
operator!=(MatRange const& a, MatRange const& b) { return !operator==(a,b); }

bool inline
isTransposed(MatRange const& i) { return (i.rs==i.cn && i.cs==1); }

bool inline
isNormal(MatRange const& i) { return (i.rs==1 && i.cs==i.rn); }

bool inline
isContiguous(MatRange const& i) { return isNormal(i) || isTransposed(i); }

inline std::ostream&
operator<<(std::ostream& s, MatRange const& mr)
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

MatRange::iterator inline MatRange::
begin() const { return iterator(*this); }

MatRange::iterator inline MatRange::
end() const { return iterator::makeEnd(*this); }

} //namespace itensor

#endif
