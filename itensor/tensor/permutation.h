//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_PERMUTATION_H
#define __ITENSOR_PERMUTATION_H

#include "itensor/util/infarray.h"

namespace itensor {

//
// Tell where each index will go, 
// if(p.dest(2) == 1) then 2 -> 1, etc.
//
struct Permutation
    {
    using storage = InfArray<long,11ul>;
    using size_type = long;
    using iterator = storage::iterator;
    using const_iterator = storage::const_iterator;
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

    storage const& 
    store() const { return store_; }

    const_iterator
    begin() const { return store_.begin(); }

    const_iterator
    end() const { return store_.end(); }

    };

inline Permutation::
Permutation(size_type size) 
  : store_(size)
    { 
    for(size_type n = 0; n < size; ++n)
        store_[n] = n;
    }


Permutation inline
inverse(Permutation const& P)
    {
    auto inv = Permutation(P.size());
    for(decltype(P.size()) n = 0; n < P.size(); ++n) 
        inv.setFromTo(P.dest(n),n);
    return inv;
    }

bool inline
isTrivial(Permutation const& P)
    {
    for(decltype(P.size()) n = 0; n < P.size(); ++n) 
        {
        if(P.dest(n) != n) return false;
        }
    return true;
    }

template<typename Container>
void
permute(Permutation const& P,
        Container const& from,
        Container & to)
    {
    for(decltype(from.size()) i = 0; i < from.size(); ++i)
        {
        to[P.dest(i)] = from[i];
        }
    }


template <typename Set1,
          typename Set2>
void
calcPerm(Set1 const& s1,
         Set2 const& s2,
         Permutation & P)
    {
    for(decltype(s2.size()) i2 = 0; i2 < s2.size(); ++i2)
        {
        auto& v2 = s2[i2];
        decltype(s1.size()) i1 = 0;
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

template <typename Set1,
          typename Set2>
Permutation
calcPerm(Set1 const& s1,
         Set2 const& s2)
    {
    auto P = Permutation(s1.size());
    calcPerm(s1,s2,P);
    return P;
    }

inline std::ostream& 
operator<<(std::ostream & s, Permutation const& P)
    {
    for(decltype(P.size()) i = 0; i < P.size(); ++i) 
        s << "(" << i << "," << P.dest(i) << ")";
    return s;
    }

} //namespace itensor

#endif
