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
#ifndef __ITENSOR_TENITER_H_
#define __ITENSOR_TENITER_H_

#include "itensor/tensor/rangeiter.h"

namespace itensor {


template<typename range_type>
auto
rangeBegin(range_type const& r) -> decltype(r.begin())
    {
    return r.begin();
    }

template<typename range_type>
auto
rangeEnd(range_type const& r) -> decltype(r.end())
    {
    return r.end();
    }

template<class Ptr, class RangeT>
class TenIter
    { 
    public:
    using value_type = typename std::iterator_traits<Ptr>::value_type;
    using reference = typename std::iterator_traits<Ptr>::reference;
    using difference_type = typename std::iterator_traits<Ptr>::difference_type;
    using pointer = typename std::iterator_traits<Ptr>::pointer;
    using iterator_category = std::forward_iterator_tag;
    using range_type = RangeT;
    using range_iter = decltype(rangeBegin(std::declval<range_type>()));
    using storage_type = DataRange<typename std::remove_pointer<Ptr>::type>;
    private:
    storage_type d_;
    range_iter it_;
    public: 

    TenIter() { }

    TenIter(storage_type d, range_type const& r) 
      : d_(d), 
        it_(rangeBegin(r))
        { }  

    reference 
    operator*() 
        { 
        return d_[it_.offset()];
        }

    TenIter& 
    operator++() { increment(); return *this; } 

    TenIter 
    operator++(int) { auto ct = *this; ct.increment(); return ct; } 

    //pointer
    //data() const { return p_; }

    range_type const&
    range() const { return it_.range(); }

    bool
    notDone() const { return it_.notDone(); }

    bool
    operator!=(TenIter const& other) const { return it_!=other.it_; }

    bool
    operator==(TenIter const& other) const { return it_==other.it_; }

    private:

    void
    increment() { ++it_; }

    public:
    //For developer use only; for making end iterator
    TenIter static
    makeEnd(range_type const& r) 
        {
        TenIter end;
        end.it_ = rangeEnd(r);
        return end;
        }
    }; 

} //namespace itensor

#endif

