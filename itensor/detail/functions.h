//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_FUNCTIONS_H
#define __ITENSOR_FUNCTIONS_H

namespace itensor {
namespace detail {

template <typename Set1,
          typename Set2,
          typename Perm>
void
calc_permutation(const Set1& s1,
                 const Set2& s2,
                 Perm& p)
    {
    for(size_t i2 = 0; i2 < s2.size(); ++i2)
        {
        const auto& v2 = s2[i2];
        bool found = false;
        for(size_t i1 = 0; i1 < s1.size(); ++i1)
            {
            if(v2 == s1[i1])
                {
                p[i1] = i2;
                found = true;
                break;
                }
            }

        if(!found)
            {
            throw ITError("sets are not permutations of each other");
            }
        }
    }


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
        const auto& v2 = *it;
        bool found = false;
        for(size_t i1 = 0; i1 < s1.size(); ++i1)
            {
            if(v2 == s1[i1])
                {
                r[i1] = m(v2);
                found = true;
                break;
                }
            }

        if(!found)
            {
            throw ITError("sets are not permutations of each other");
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
Real inline
quickran()
    {
    static int seed = (std::time(NULL) + getpid());
    int im = 134456;
    int ia = 8121;
    int ic = 28411;
    Real scale = 1.0 / im;
    seed = (seed*ia+ic)%im;
    return Real(seed) * scale;
    }


}; //namespace detail
}; //namespace itensor

#endif

