//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "vector.h"

namespace itensor {

void vecref::
operator=(const vecref& other)
    {
    store_ = other.store_;
    cstore_ = other.cstore_;
    strd_ = other.strd_;
    size_ = other.size_;
    }


std::ostream&
operator<<(std::ostream& s, const vecref& v)
    {
    for(auto j = 1l; j <= v.size(); ++j)
        {
        s << v(j) << " ";
        }
    return s;
    }

Real
norm(const vecref& v)
    {
    Real nrm = 0;
    for(auto& el : v) nrm += el*el;
    return std::sqrt(nrm);
    }

}; //namespace itensor
