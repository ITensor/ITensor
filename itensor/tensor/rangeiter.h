//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_RANGEITER_H_
#define __ITENSOR_RANGEITER_H_

#include <iterator> 
#include "range.h"

namespace itensor {

template<typename range_type_>
class RangeIter
    { 
    public:
    using range_type = range_type_;
    using size_type = typename range_type::size_type;
    using value_type = size_type;
    using reference = size_type &;
    //using difference_type = typename std::iterator_traits<value_type>::difference_type;
    //using pointer = typename std::iterator_traits<value_type>::pointer;
    using iterator_category = std::forward_iterator_tag;
    using ind_type = InfArray<size_type,8ul>;
    using const_ind_iterator = typename ind_type::const_iterator;
    private:
    const range_type *prange_ = nullptr; 
    value_type val_ = 0;
    ind_type ind_;
    public: 

    RangeIter()
      : prange_(nullptr),
        val_(0)
        { }

    RangeIter(range_type const& R) 
      : prange_(&R),
        val_(0),
        ind_(R.r(),0)
        { }

    RangeIter(RangeIter && o) { operator=(std::move(o)); }

    RangeIter&
    operator=(RangeIter && o)
        {
        prange_ = o.prange_;
        val_ = o.val_;
        ind_ = std::move(o.ind_);
        return *this;
        }

    RangeIter(RangeIter const& o) = delete;

    RangeIter&
    operator=(RangeIter const& o) = delete;

    RangeIter const&
    operator*() const { return *this; }  

    size_type
    operator[](size_type n) const { return ind_[n]; }

    size_type
    size() const { return ind_.size(); }

    value_type
    val() const { return val_; }

    range_type const&
    range() const { return *prange_; }

    RangeIter& 
    operator++() { increment(); return *this; } 

    RangeIter 
    operator++(int) { auto ct = *this; ct.increment(); return ct; } 

    ind_type
    ind() const { return ind_; }

    bool
    operator==(RangeIter const& o) const 
        { 
#ifdef DEBUG
        if(prange_ != o.prange_) 
            Error("Can't compare RangeIter created from different range objects");
#endif
        return val_==o.val_;
        }

    bool
    operator!=(RangeIter const& o) const { return !operator==(o); }

    const_ind_iterator
    begin() const { return ind_.begin(); }

    const_ind_iterator
    end() const { return ind_.end(); }

    private:

    void
    increment()
        {
#ifdef DEBUG
        if(range().r() == 0) Error("Can't increment RangeIter made from rank 0 range");
#endif
        auto r = range().r();
        ind_[0] += 1;
        val_ += range().stride(0);
        if(ind_[0] == range().extent(0))
            {
            for(decltype(r) n = 1; n < r; ++n)
                {
                ind_[n-1] = 0;
                val_ -= range().extent(n-1)*range().stride(n-1);
                ind_[n] += 1;
                val_ += range().stride(n);
                if(ind_[n] < range().extent(n)) return;
                }
            //will only reach this line when totally done
            val_ = std::numeric_limits<value_type>::max();
            }
        }
    public:
    //For developer use only; for making end iterator
    RangeIter static
    makeEnd(range_type const& R) 
        {
        RangeIter end;
        end.val_ = std::numeric_limits<value_type>::max();
        end.prange_ = &R; 
        return end;
        }
    }; 

template<typename R>
std::ostream&
operator<<(std::ostream & s,
           RangeIter<R> const& it)
    {
    auto r = it.range().r();
    s << format("%*d",3,it.val()) << " (";
    for(decltype(r) j = 0; j < r; ++j)
        {
        s << it[j];
        if(1+j != r) s << ",";
        }
    s << ")";
    return s;
    }

} //namespace itensor

#endif

