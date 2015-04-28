//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#include "vec.h"
#include "lapack_wrap.h"
#include <limits>
#include "detail/algs.h"

namespace itensor {

//
//
//  \       /  /----   /-----  .----.   /----  /----
//   \     /   |       |       |    |   |      |    
//    \   /    |---    |       |---'.   |---   |--- 
//     \ /     |       |       |    |   |      |    
//      ^      \----   \-----  |    |   \----  |    
//
//

template<typename Func, typename Iter>
void
apply(VecRef& v,
      Iter it,
      const Func& f)
    {
    for(auto& el : v) 
        {
        f(el,*it);
        ++it;
        }
    }

VecRef
operator&=(VecRef a, VecRefc b)
    {
#ifdef DEBUG
    if(b.size() != a.size()) throw std::runtime_error("mismatched sizes in VecRef operator&=");
#endif
    auto assign = [](Real& x, Real y) { x = y; };
    if(b.contiguous()) apply(a,b.data(),assign);
    else               apply(a,b.cbegin(),assign);
    return a;
    }

VecRef 
operator*=(VecRef a, Real fac)
    {
    if(a.contiguous())
        {
#ifdef DEBUG
        if(a.size() > std::numeric_limits<LAPACK_INT>::max()) 
            throw std::runtime_error("VecRef overflow of size beyond LAPACK_INT range");
#endif
        dscal_wrapper(a.size(),fac,a.data());
        return a;
        }
    for(auto& el : a) el *= fac;
    return a;
    }

VecRef 
operator/=(VecRef a, Real fac)
    {
    if(fac == 0) throw std::runtime_error("VecRef /=: divide by zero");
    return operator*=(a,1./fac);
    }

void
call_daxpy(VecRef& A, const VecRefc& B, Real alpha_)
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

VecRef
operator+=(VecRef a, VecRefc b)
    {
#ifdef DEBUG
    if(a.size()!=b.size()) throw std::runtime_error("VecRef +=: mismatched sizes");
#endif
    auto pluseq = [](Real& x, Real y) { x += y; };
    if(b.contiguous())
        {
        if(a.contiguous()) call_daxpy(a,b,+1);
        else               apply(a,b.data(),pluseq);
        }
    else
        {
        apply(a,b.cbegin(),pluseq);
        }
    return a;
    }

VecRef
operator-=(VecRef a, VecRefc b)
    {
#ifdef DEBUG
    if(a.size()!=b.size()) throw std::runtime_error("VecRef +=: mismatched sizes");
#endif
    auto minuseq = [](Real& x, Real y) { x -= y; };
    if(b.contiguous())
        {
        if(a.contiguous()) call_daxpy(a,b,-1);
        else               apply(a,b.data(),minuseq);
        }
    else
        {
        apply(a,b.cbegin(),minuseq);
        }
    return a;
    }

std::ostream&
operator<<(std::ostream& s, VecRefc v)
    {
    for(const auto& el : v)
        {
        s << el << " ";
        }
    return s;
    }

Real
norm(VecRefc v)
    {
    Real nrm = 0;
    for(auto& el : v) nrm += el*el;
    return std::sqrt(nrm);
    }

Real
operator*(VecRefc a, VecRefc b)
    {
#ifdef DEBUG
    if(a.size() != b.size()) throw std::runtime_error("VecRef dot product: mismatched sizes");
    if(a.size() > std::numeric_limits<LAPACK_INT>::max()) 
        throw std::runtime_error("VecRef dot product: overflow of size beyond LAPACK_INT range");
#endif
    return ddot_wrapper(a.size(),a.data(),a.stride(),b.data(),b.stride());
    }

void
randomize(VecRef v)
    {
    for(auto& el : v) el = detail::quickran();
    }

//
//
//  \       /   /----   /-----
//   \     /    |       |
//    \   /     |---    |
//     \ /      |       |
//      ^       \----   \-----
//
//


Vec
randomVec(long size)
    {
    Vec v(size);
    randomize(v);
    return v;
    }


//bool
//overlaps(const Real* b1, const Real* e1,
//         const Real* b2, const Real* e2)
//    {
//    using ptr_less = std::less<const Real*>;
//    return (ptr_less(b1,b2) && ptr_less(b2,e1))
//        || (ptr_less(b1,e2) && ptr_less(e2,e1));
//    }


}; //namespace itensor
