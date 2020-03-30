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
#ifndef __ITENSOR_ALGS_H
#define __ITENSOR_ALGS_H

#include <sys/types.h>
#include <unistd.h>
#include <algorithm>
#include <ctime>
#include <stdexcept>
#include <random>
#include "itensor/types.h"

namespace itensor {
namespace detail {



template <typename Set1,
          typename Set2Iter,
          typename RType,
          typename Map>
void
permute_map(const Set1& s1,
            const Set2Iter& s2begin,
            const Set2Iter& s2end,
            RType& r,
            Map&& m)
    {
    for(auto it = s2begin; it != s2end; ++it)
        {
        auto& v2 = *it;
        bool found = false;
        for(size_t i1 = 0; i1 < s1.size(); ++i1)
            if(v2 == s1[i1])
                {
                r[i1] = m(v2);
                found = true;
                break;
                }

        if(!found)
            {
            throw std::runtime_error("sets are not permutations of each other");
            }
        }
    }

template <typename Set1,
          typename Set2,
          typename RType,
          typename Map>
void
permute_map(const Set1& s1,
            const Set2& s2,
            RType& r,
            Map&& m)
    {
    permute_map(s1,std::begin(s2),std::end(s2),r,std::forward<Map>(m));
    }

template <typename Container, typename Item>
bool
contains(const Container& C,
         const Item& I)
    {
    for(const auto& c : C) 
        {
        if(I == c) return true;
        }
    return false;
    }

//Simple linear congruential random number generator
inline long&
seed_quickran(int newseed)
    {
    static long seed = (std::time(NULL) + getpid());
    if(newseed != 0) seed = newseed;
    return seed;
    }

double inline
quickran() 
    {
    auto res = 0.0;
    while(res == 0.0)
        {
        long& seed = seed_quickran(0);
        long im = 134456;
        long ia = 8121;
        long ic = 28411;
        double scale = 1.0 / im;
        seed = (seed*ia+ic)%im;
        res = std::fabs(double(seed)) * scale;
        }
    return res;
    }

Cplx inline
quickranCplx() { return Cplx(detail::quickran(),detail::quickran()); }

template<typename Container>
using const_correct_ptr = 
    typename std::conditional<std::is_const<Container>::value,
              typename Container::const_pointer,
              typename Container::pointer>::type;

template<typename Container, typename T, typename Compare>
const_correct_ptr<Container>
binaryFind(Container& c, const T& val, Compare less)
    {
    auto range = std::equal_range(std::begin(c),std::end(c),val,less);
    if(range.first != range.second)
        return &(*range.first);
    else
        return nullptr;
    }

template<typename Container, typename T>
const_correct_ptr<Container>
binaryFind(Container& c, const T& val)
    {
    auto range = std::equal_range(std::begin(c),std::end(c),val);
    if(range.first != range.second)
        return &(*range.first);
    else
        return nullptr;
    }


} //namespace detail
} //namespace itensor

#endif

