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
#ifndef __ITENSOR_ITERATE_H_
#define __ITENSOR_ITERATE_H_

#include <type_traits>
#include "itensor/util/itertools.h"

namespace itensor {

using triqs::utility::enumerate;
using triqs::utility::zip;
using triqs::utility::product;

namespace detail {

template <typename size_type>
class RangeHelper
    {
    size_type curr_;
    size_type end_;
    public:

    constexpr
    RangeHelper(size_type b,
                size_type e)
      : curr_(b),
        end_(e)
        { }

    size_type const&
    operator*() const { return curr_; }

    RangeHelper& 
    operator++() 
        { 
        ++curr_; 
        return *this;
        }

    bool
    operator!=(RangeHelper const& other) const
        {
        return curr_ < other.curr_;
        }

    RangeHelper 
    begin() const { return RangeHelper(curr_,end_); }

    RangeHelper 
    end() const { return RangeHelper(end_,end_); }
    };

} //namespace detail

//
//  No way to get the first index, so we can only support current.
//
template <typename T> inline
T
current(const detail::RangeHelper<T>& rh)
{
    return *rh;
}
template <typename T> inline
T
last(const detail::RangeHelper<T>& rh)
{
    return *rh.end();
}
template <typename T> inline
T
length(const detail::RangeHelper<T>& rh)
{
    return last(rh)-current(rh);
}

template <typename T,
          class=typename std::enable_if<std::is_integral<T>::value>::type>
auto constexpr
range(T end) 
    -> detail::RangeHelper<T>
    {
    return detail::RangeHelper<T>(0,end);
    }

template <typename ST, typename T,
          class=typename std::enable_if<
                         std::is_integral<ST>::value
                      && std::is_integral<T>::value>::type>
auto constexpr
range(ST start, T end) 
    -> detail::RangeHelper<T>
    {
    return detail::RangeHelper<T>(T(start),end);
    }

template <typename T,
          class=typename std::enable_if<std::is_integral<T>::value>::type>
auto constexpr
range1(T end) -> detail::RangeHelper<T>
    {
    return detail::RangeHelper<T>(1,1+end);
    }

template <typename ST, typename T,
          class=typename std::enable_if<
                         std::is_integral<ST>::value
                      && std::is_integral<T>::value>::type>
auto constexpr
range1(ST start, T end) -> detail::RangeHelper<T>
    {
    return detail::RangeHelper<T>(start,1+end);
    }

template <typename C> constexpr
auto
range(C const& container) 
    -> detail::RangeHelper<decltype(container.size())>
    {
    using size_type = decltype(container.size());
    return detail::RangeHelper<size_type>(0,container.size());
    }

template <typename C> constexpr
auto
range1(C const& container) 
    -> detail::RangeHelper<decltype(container.size())>
    {
    using size_type = decltype(container.size());
    return detail::RangeHelper<size_type>(1,1+container.size());
    }

}

#endif
