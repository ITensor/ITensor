//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#include "itensor/matrix/vec.h"
#include <limits>
#include "itensor/detail/algs.h"

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
apply(VectorRef& v,
      Iter it,
      Func const& f)
    {
    for(auto& el : v) 
        {
        f(el,*it);
        ++it;
        }
    }

VectorRef
operator&=(VectorRef a, VectorRefc b)
    {
#ifdef DEBUG
    if(b.size() != a.size()) throw std::runtime_error("mismatched sizes in VectorRef operator&=");
#endif
    auto assign = [](Real& x, Real y) { x = y; };
    if(isContiguous(b)) apply(a,b.data(),assign);
    else                apply(a,b.cbegin(),assign);
    return a;
    }

VectorRef 
operator*=(VectorRef a, Real fac)
    {
    if(isContiguous(a))
        {
#ifdef DEBUG
        if(a.size() > std::numeric_limits<LAPACK_INT>::max()) 
            throw std::runtime_error("VectorRef overflow of size beyond LAPACK_INT range");
#endif
        dscal_wrapper(a.size(),fac,a.data());
        return a;
        }
    for(auto& el : a) el *= fac;
    return a;
    }

VectorRef 
operator/=(VectorRef a, Real fac)
    {
    if(fac == 0) throw std::runtime_error("VectorRef /=: divide by zero");
    return operator*=(a,1./fac);
    }

void
call_daxpy(VectorRef& A, const VectorRefc& B, Real alpha_)
    {
    LAPACK_REAL alpha = alpha_;
    LAPACK_INT inc = 1;
    LAPACK_INT size = A.size();
#ifdef DEBUG
    if(A.size() > std::numeric_limits<LAPACK_INT>::max()) 
        throw std::runtime_error("overflow of size beyond LAPACK_INT range");
#endif
    daxpy_wrapper(size,alpha,B.data(),inc,A.data(),inc);
    }

VectorRef
operator+=(VectorRef a, VectorRefc b)
    {
#ifdef DEBUG
    if(a.size()!=b.size()) throw std::runtime_error("VectorRef +=: mismatched sizes");
#endif
    auto pluseq = [](Real& x, Real y) { x += y; };
    if(isContiguous(b))
        {
        if(isContiguous(a)) call_daxpy(a,b,+1);
        else                apply(a,b.data(),pluseq);
        }
    else
        {
        apply(a,b.cbegin(),pluseq);
        }
    return a;
    }

VectorRef
operator-=(VectorRef a, VectorRefc b)
    {
#ifdef DEBUG
    if(a.size()!=b.size()) throw std::runtime_error("VectorRef +=: mismatched sizes");
#endif
    auto minuseq = [](Real& x, Real y) { x -= y; };
    if(isContiguous(b))
        {
        if(isContiguous(a)) call_daxpy(a,b,-1);
        else                apply(a,b.data(),minuseq);
        }
    else
        {
        apply(a,b.cbegin(),minuseq);
        }
    return a;
    }

template<>
std::ostream&
operator<<(std::ostream& s, VectorRefc const& v)
    {
    for(auto& el : v) s << el << " ";
    return s;
    }

template<>
Real
norm(VectorRefc const& v)
    {
    return dnrm2_wrapper(v.size(),v.data(),stride(v)); 
    }

Real
operator*(VectorRefc a, VectorRefc b)
    {
#ifdef DEBUG
    if(a.size() != b.size()) throw std::runtime_error("VectorRef dot product: mismatched sizes");
    if(a.size() > std::numeric_limits<LAPACK_INT>::max()) 
        throw std::runtime_error("VectorRef dot product: overflow of size beyond LAPACK_INT range");
#endif
    return ddot_wrapper(a.size(),a.data(),stride(a),b.data(),stride(b));
    }

void
randomize(VectorRef v)
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


Vector
randomVec(long size)
    {
    Vector v(size);
    randomize(v);
    return v;
    }

Real
sumels(VectorRefc v)
    {
    Real tot = 0;
    for(auto& el : v) tot += el;
    return tot;
    }

void
resize(Vector & v,
       size_t newsize)
    {
    v.resize(VecRange(newsize));
    }


//bool
//overlaps(const Real* b1, const Real* e1,
//         const Real* b2, const Real* e2)
//    {
//    using ptr_less = std::less<const Real*>;
//    return (ptr_less(b1,b2) && ptr_less(b2,e1))
//        || (ptr_less(b1,e2) && ptr_less(e2,e1));
//    }


} //namespace itensor
