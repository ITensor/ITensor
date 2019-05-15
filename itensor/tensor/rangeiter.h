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
#ifndef __ITENSOR_RANGEITER_H_
#define __ITENSOR_RANGEITER_H_

#include <iterator> 
#include "itensor/tensor/types.h"

namespace itensor {

template<typename range_type_>
class RangeIter
    { 
    public:
    using range_type = range_type_;
    using size_type = typename range_type::size_type;
    using offset_type = size_type;
    using reference = size_type &;
    using iterator_category = std::forward_iterator_tag;
    using ind_type = InfArray<size_type,11ul>;
    using const_ind_iterator = typename ind_type::const_iterator;
    private:
    range_type const* prange_ = nullptr; 
    offset_type off_ = 0;
    ind_type ind_;
    public: 

    RangeIter()
      : prange_(nullptr),
        off_(0)
        { }

    RangeIter(range_type const& R) 
      : prange_(&R),
        off_(0),
        ind_(R.order(),R.start())
        { }


    RangeIter const&
    operator*() const { return *this; }  

    size_type const&
    operator[](size_type n) const { return ind_[n]; }

    size_type
    size() const { return ind_.size(); }

    offset_type
    offset() const { return off_; }

    ind_type
    index() const { return ind_; }

    range_type const&
    range() const { return *prange_; }

    bool
    notDone() const
        {
        return off_ != std::numeric_limits<offset_type>::max();
        }

    RangeIter& 
    operator++() { increment(); return *this; } 

    RangeIter 
    operator++(int) { auto ct = *this; ct.increment(); return ct; } 


    bool
    operator==(RangeIter const& o) const 
        { 
#ifdef DEBUG
        if(prange_ != o.prange_) 
            Error("Can't compare RangeIter created from different range objects");
#endif
        return off_==o.off_;
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
        using rextent = decltype(range().extent(0));
#ifdef DEBUG
        if(range().order() == 0) Error("Can't increment RangeIter made from order 0 range");
#endif
        auto r = range().order();
        ind_[0] += 1;
        off_ += range().stride(0);
        if(rextent(ind_[0]-range().start()) == range().extent(0))
            {
            for(decltype(r) n = 1; n < r; ++n)
                {
                ind_[n-1] = range().start();
                off_ -= range().extent(n-1)*range().stride(n-1);
                ind_[n] += 1;
                off_ += range().stride(n);
                if(rextent(ind_[n]-range().start()) < range().extent(n)) return;
                }
            //will only reach this line when totally done
            off_ = std::numeric_limits<offset_type>::max();
            }
        }
    public:
    //For developer use only; for making end iterator
    RangeIter static
    makeEnd(range_type const& R) 
        {
        RangeIter end;
        end.off_ = std::numeric_limits<offset_type>::max();
        end.prange_ = &R; 
        return end;
        }
    }; 

template<typename R>
std::ostream&
operator<<(std::ostream & s,
           RangeIter<R> const& it)
    {
    auto r = it.range().order();
    s << format("%*d",3,it.offset()) << " (";
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

