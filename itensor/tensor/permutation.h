//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_PERMUTATION_H
#define __ITENSOR_PERMUTATION_H

#include "infarray.h"

namespace itensor {

//
// Tell where each index will go, 
// if(p.dest(2) == 1) then 2 -> 1, etc.
//
struct Permutation
    {
    using storage = InfArray<long,8>;
    using size_type = long;
    private: 
    storage store_;
    public:

    Permutation() { }

    Permutation(size_type size);

    size_type
    size() const { return store_.size(); }

    explicit operator bool() const { return !store_.empty(); }

    void 
    setFromTo(size_type from, 
              size_type to)
        {
        store_[from] = to; 
        }

    long 
    dest(size_type j) const { return store_[j]; }

    long&
    operator[](size_type j) { return store_[j]; }

    long
    operator[](size_type j) const { return store_[j]; }

    const storage& 
    store() const { return store_; }

    };

inline Permutation::
Permutation(size_type size) 
    : 
    store_(size)
    { 
    for(size_type n = 0; n < size; ++n)
        store_[n] = n;
    }


Permutation inline
inverse(const Permutation& P)
    {
    Permutation inv(P.size());
    for(Permutation::size_type n = 0; n < P.size(); ++n) 
        inv.setFromTo(P.dest(n),n);
    return inv;
    }

bool inline
isTrivial(const Permutation& P)
    {
    for(Permutation::size_type n = 0; n < P.size(); ++n) 
        {
        if(P.dest(n) != n) return false;
        }
    return true;
    }

template<typename Container>
void
permute(const Permutation& P,
        const Container& from,
        Container& to)
    {
    using size_type = typename Container::size_type;
    for(size_type i = 0; i < from.size(); ++i)
        {
        to[P.dest(i)] = from[i];
        }
    }

template <typename Set1,
          typename Set2>
void
calc_permutation(const Set1& s1,
                 const Set2& s2,
                 Permutation& P)
    {
    using size_type1 = decltype(s1.size());
    using size_type2 = decltype(s2.size());
    for(size_type2 i2 = 0; i2 < s2.size(); ++i2)
        {
        auto& v2 = s2[i2];
        size_type1 i1 = 0;
        for(; i1 < s1.size(); ++i1)
            {
            if(v2 == s1[i1])
                {
                P[i1] = i2;
                break;
                }
            }
        if(i1 == s1.size())
            {
            throw std::runtime_error("sets are not permutations of each other");
            }
        }
    }

inline std::ostream& 
operator<<(std::ostream& s, const Permutation& P)
    {
    for(Permutation::size_type i = 0; i < P.size(); ++i) 
        s << "(" << i << "," << P.dest(i) << ")";
    return s;
    }

} //namespace itensor

#endif
