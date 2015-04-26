//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#include "mat.h"
#include "lapack_wrap.h"
#include <limits>
#include "detail/algs.h"

namespace itensor {

//
//
//  .      .    __    _______   ____    .____  .____
//  |\    /|   /  \      |     |    |   |      |    
//  | \  / |  |____|     |     |___/    |___   |--- 
//  |  \/  |  |    |     |     |   \    |      |    
//  |      |  |    |     |     |    |   |____  |    
//
//

template<typename Func, typename Iter>
void
apply(MatRef& v,
      Iter it,
      const Func& f)
    {
    for(auto& el : v) 
        {
        f(el,*it);
        ++it;
        }
    }

MatRef& 
operator&=(MatRef& a, CMatRef b)
    {
#ifdef DEBUG
    if(!(b.Nrows()==a.Nrows() && b.Ncols()==a.Ncols())) 
        throw std::runtime_error("mismatched sizes in VecRef operator&=");
#endif
    auto assign = [](Real& x, Real y) { x = y; };
    if(b.contiguous()) apply(a,b.data(),assign);
    else               apply(a,b.cbegin(),assign);
    return a;
    }

MatRef& 
operator*=(MatRef& a, Real fac)
    {
    if(a.contiguous())
        {
#ifdef DEBUG
        if(a.size() > std::numeric_limits<LAPACK_INT>::max()) 
            throw std::runtime_error("MatRef overflow of size beyond LAPACK_INT range");
#endif
        dscal_wrapper(a.size(),fac,a.data());
        }
    else
        {
        for(auto& el : a) el *= fac;
        }
    return a;
    }

MatRef& 
operator/=(MatRef& a, Real fac)
    {
    if(fac == 0) throw std::runtime_error("MatRef /=: divide by zero");
    return operator*=(a,1./fac);
    }

void
call_daxpy(MatRef& A, const CMatRef& B, Real alpha_)
    {
    LAPACK_REAL alpha = alpha_;
    LAPACK_INT inc = 1;
    LAPACK_INT size = A.size();
#ifdef DEBUG
    if(A.size() > std::numeric_limits<LAPACK_INT>::max()) 
        throw std::runtime_error("overflow of size beyond LAPACK_INT range");
#endif
    daxpy_wrapper(&size,&alpha,B.data(),&inc,A.data(),&inc);
    }

MatRef&
operator+=(MatRef& a, CMatRef b)
    {
#ifdef DEBUG
    if(!(a.Ncols()==b.Ncols() && a.Nrows()==b.Nrows())) 
        throw std::runtime_error("MatRef +=: mismatched sizes");
#endif
    if(b.ind()==a.ind() && b.contiguous())
        {
        call_daxpy(a,b,+1);
        }
    else
        {
        auto pluseq = [](Real& x, Real y) { x += y; };
        apply(a,b.cbegin(),pluseq);
        }
    return a;
    }

MatRef&
operator-=(MatRef& a, CMatRef b)
    {
#ifdef DEBUG
    if(!(a.Ncols()==b.Ncols() && a.Nrows()==b.Nrows())) 
        throw std::runtime_error("MatRef +=: mismatched sizes");
#endif
    if(b.ind()==a.ind() && b.contiguous())
        {
        call_daxpy(a,b,-1);
        }
    else
        {
        auto minuseq = [](Real& x, Real y) { x -= y; };
        apply(a,b.cbegin(),minuseq);
        }
    return a;
    }

Real
norm(CMatRef v)
    {
    Real nrm = 0;
    for(auto& el : v) nrm += el*el;
    return std::sqrt(nrm);
    }


std::ostream&
operator<<(std::ostream& s, CMatRef M)
    {
    for(long r = 1; r <= M.Nrows(); ++r)
        {
        s << "|";
        for(long c = 1; c <= M.Ncols(); ++c)
            {
            s << M(r,c);
            s << (c == M.Ncols() ? "|" : " ");
            }
        if(r < M.Nrows()) s << "\n";
        }
    return s;
    }


}; //namespace itensor
