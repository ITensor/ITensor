//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "vector.h"
#include "lapack_wrap.h"
#include <limits>

namespace itensor {

void vecref::
operator=(const vecref& other)
    {
    store_ = other.store_;
    cstore_ = other.cstore_;
    strd_ = other.strd_;
    size_ = other.size_;
    }

void vecref::
operator*=(Real fac)
    {
#ifdef DEBUG
    if(readOnly()) throw std::runtime_error("vecref *=: read only");
#endif
    if(contiguous())
        {
#ifdef DEBUG
        if(size() > std::numeric_limits<LAPACK_INT>::max()) 
            throw std::runtime_error("matrixref *=: overflow of size beyond LAPACK_INT range");
#endif
        dscal_wrapper(size(),fac,store());
        }
    else
        {
        for(auto& el : *this) el *= fac;
        }
    }
void vecref::
operator/=(Real fac)
    {
    if(fac == 0) throw std::runtime_error("vecref /=: divide by zero");
    operator*=(1./fac);
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
