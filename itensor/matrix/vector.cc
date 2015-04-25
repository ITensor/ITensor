//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "vector.h"
#include "lapack_wrap.h"
#include <limits>
#include "detail/algs.h"

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

void
call_daxpy(vecref& A, const vecref& B, Real alpha_)
    {
    LAPACK_REAL alpha = alpha_;
    LAPACK_INT inc = 1;
    LAPACK_INT size = A.size();
#ifdef DEBUG
    if(A.size() > std::numeric_limits<LAPACK_INT>::max()) throw std::runtime_error("overflow of size beyond LAPACK_INT range");
#endif
    daxpy_wrapper(&size,&alpha,B.cstore(),&inc,A.store(),&inc);
    }

void vecref::
operator+=(const vecref& other)
    {
#ifdef DEBUG
    if(size()!=other.size()) throw std::runtime_error("vecref+=: mismatched sizes");
#endif
    if(contiguous() && other.contiguous())
        {
        call_daxpy(*this,other,+1);
        }
    else
        {
        auto o = other.begin();
        for(auto& el : *this) 
            {
            el += *o;
            ++o;
            }
        }
    }
void vecref::
operator-=(const vecref& other)
    {
#ifdef DEBUG
    if(size()!=other.size()) throw std::runtime_error("vecref+=: mismatched sizes");
#endif
    if(contiguous() && other.contiguous())
        {
        call_daxpy(*this,other,-1);
        }
    else
        {
        auto o = other.begin();
        for(auto& el : *this) 
            {
            el -= *o;
            ++o;
            }
        }
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

vec
randomVec(long size)
    {
    vec res(size);
    for(auto& el : res) el = detail::quickran();
    return res;
    }

}; //namespace itensor
