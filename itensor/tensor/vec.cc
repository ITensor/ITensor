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
#include <limits>
#include "itensor/detail/algs.h"
#include "itensor/tensor/vec.h"
#include "itensor/tensor/lapack_wrap.h"

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

template<typename T, typename Func, typename Iter>
void
apply(VecRef<T> const& v,
      Iter it,
      Func const& f)
    {
    for(auto& el : v) 
        {
        f(el,*it);
        ++it;
        }
    }

template<typename T1, typename T2>
void
refAssign(VecRef<T1> const& a, VecRefc<T2> const& b)
    {
#ifdef DEBUG
    if(b.size() != a.size()) throw std::runtime_error("mismatched sizes in VectorRef operator&=");
#endif
    auto assign = [](T1& x, T2 y) { x = y; };
    if(isContiguous(b)) apply(a,b.data(),assign);
    else                apply(a,b.cbegin(),assign);
    }

void
operator&=(VectorRef a, VectorRefc const& b)
    {
    refAssign(a,b);
    }
void
operator&=(CVectorRef a, CVectorRefc const& b)
    {
    refAssign(a,b);
    }
void
operator&=(CVectorRef a, VectorRefc const& b)
    {
    refAssign(a,b);
    }

template<typename V>
void 
multReal(VecRef<V> const& v, Real fac)
    {
    if(isContiguous(v))
        {
#ifdef DEBUG
        if(v.size() > std::numeric_limits<unsigned long>::max()) 
            throw std::runtime_error("VectorRef overflow of size beyond long unsigned int range");
#endif
        auto d = realData(v);
        dscal_wrapper(d.size(),fac,d.data());
        }
    else
        {
        for(auto& el : v) el *= fac;
        }
    }

VectorRef
operator*=(VectorRef v, Real fac)
    {
    multReal(v,fac);
    return v;
    }

CVectorRef
operator*=(CVectorRef v, Real fac)
    {
    multReal(v,fac);
    return v;
    }

CVectorRef
operator*=(CVectorRef v, Cplx fac)
    {
    for(auto& el : v) el *= fac;
    return v;
    }

template<typename T>
void
divReal(VecRef<T> V, Real fac)
    {
    if(fac == 0) throw std::runtime_error("VectorRef /=: divide by zero");
    if(isContiguous(V))
        {
        auto d = realData(V);
        auto vend = d.data()+d.size();
        for(auto i = d.data(); i != vend; ++i)
            {
            *i /= fac;
            }
        }
    else
        {
        for(auto& el : V) el /= fac;
        }
    }

VectorRef
operator/=(VectorRef v, Real fac)
    {
    divReal(v,fac);
    return v;
    }

CVectorRef
operator/=(CVectorRef v, Real fac)
    {
    divReal(v,fac);
    return v;
    }

void
call_daxpy(VectorRef& A, const VectorRefc& B, Real alpha_)
    {
    LAPACK_REAL alpha = alpha_;
    LAPACK_INT inc = 1;
    LAPACK_INT size = A.size();
#ifdef DEBUG
    if(A.size() > std::numeric_limits<unsigned long>::max()) 
        throw std::runtime_error("overflow of size beyond long unsigned int range");
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
    if(a.size() > std::numeric_limits<unsigned long>::max()) 
        throw std::runtime_error("VectorRef dot product: overflow of size beyond long unsigned int range");
#endif
    return ddot_wrapper(a.size(),a.data(),stride(a),b.data(),stride(b));
    }

Cplx
operator*(CVectorRefc a, CVectorRefc b)
    {
#ifdef DEBUG
    if(a.size() != b.size()) throw std::runtime_error("VectorRef dot product: mismatched sizes");
    if(a.size() > std::numeric_limits<unsigned long>::max()) 
        throw std::runtime_error("VectorRef dot product: overflow of size beyond long unsigned int range");
#endif
    return zdotc_wrapper(a.size(),a.data(),stride(a),b.data(),stride(b));
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

CVector
randomCVec(long size)
    {
    CVector v(size);
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

void
resize(CVector & v,
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
