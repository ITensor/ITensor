//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "vec.h"
#include "lapack_wrap.h"
#include <limits>
#include "detail/algs.h"

namespace itensor {

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


VecRef& 
operator&=(VecRef& a, CVecRef b)
    {
#ifdef DEBUG
    if(b.size() != a.size()) throw std::runtime_error("mismatched sizes in VecRef operator&=");
#endif
    auto assign = [](Real& x, Real y) { x = y; };
    if(b.contiguous()) apply(a,b.data(),assign);
    else               apply(a,b.cbegin(),assign);
    return a;
    }

//void vecref::
//operator*=(Real fac)
//    {
//#ifdef DEBUG
//    if(readOnly()) throw std::runtime_error("vecref *=: read only");
//#endif
//    if(contiguous())
//        {
//#ifdef DEBUG
//        if(size() > std::numeric_limits<LAPACK_INT>::max()) 
//            throw std::runtime_error("matrixref *=: overflow of size beyond LAPACK_INT range");
//#endif
//        dscal_wrapper(size(),fac,store());
//        }
//    else
//        {
//        for(auto& el : *this) el *= fac;
//        }
//    }
//void vecref::
//operator/=(Real fac)
//    {
//    if(fac == 0) throw std::runtime_error("vecref /=: divide by zero");
//    operator*=(1./fac);
//    }
//
//void
//call_daxpy(vecref& A, const vecref& B, Real alpha_)
//    {
//    LAPACK_REAL alpha = alpha_;
//    LAPACK_INT inc = 1;
//    LAPACK_INT size = A.size();
//#ifdef DEBUG
//    if(A.size() > std::numeric_limits<LAPACK_INT>::max()) throw std::runtime_error("overflow of size beyond LAPACK_INT range");
//#endif
//    daxpy_wrapper(&size,&alpha,B.cstore(),&inc,A.store(),&inc);
//    }
//
//void vecref::
//operator+=(const vecref& other)
//    {
//#ifdef DEBUG
//    if(size()!=other.size()) throw std::runtime_error("vecref+=: mismatched sizes");
//#endif
//    if(contiguous() && other.contiguous())
//        {
//        call_daxpy(*this,other,+1);
//        }
//    else
//        {
//        auto o = other.begin();
//        for(auto& el : *this) 
//            {
//            el += *o;
//            ++o;
//            }
//        }
//    }
//void vecref::
//operator-=(const vecref& other)
//    {
//#ifdef DEBUG
//    if(size()!=other.size()) throw std::runtime_error("vecref+=: mismatched sizes");
//#endif
//    if(contiguous() && other.contiguous())
//        {
//        call_daxpy(*this,other,-1);
//        }
//    else
//        {
//        auto o = other.begin();
//        for(auto& el : *this) 
//            {
//            el -= *o;
//            ++o;
//            }
//        }
//    }
//
////
////  \       /   /----   /-----
////   \     /    |       |
////    \   /     |---    |
////     \ /      |       |
////      ^       \----   \-----
////
//
//void vec::
//assignFromRef(const vecref& other)
//    {
//    if(&other == this) return;
//    //Copy data from other contiguously into data_
//    data_ = storage_type(other.cbegin(),other.cend());
//    //Set up appropriate store pointers, size, stride
//    parent::operator=(vecref(data_.data(),data_.size()));
//    }
//
//void vec::
//assignFromVec(const vec& other)
//    {
//    if(&other == this) return;
//    data_ = other.data_;
//    parent::operator=(vecref(data_.data(),data_.size()));
//    }
//
//void vec::
//moveFromVec(vec&& other)
//    {
//    data_ = std::move(other.data_);
//    other.clear();
//    parent::operator=(vecref(data_.data(),data_.size()));
//    }
//
//

VecRef
randomize(VecRef v)
    {
    for(auto& el : v) el = detail::quickran();
    return v;
    }

Vec
randomVec(long size)
    {
    Vec v(size);
    randomize(v);
    return v;
    }

std::ostream&
operator<<(std::ostream& s, CVecRef v)
    {
    for(const auto& el : v)
        {
        s << el << " ";
        }
    return s;
    }

//Real
//norm(const vecref& v)
//    {
//    Real nrm = 0;
//    for(auto& el : v) nrm += el*el;
//    return std::sqrt(nrm);
//    }
//
//Real
//operator*(const vecref& A, const vecref& B)
//    {
//#ifdef DEBUG
//    if(A.size() != B.size()) throw std::runtime_error("vecref*vecref: mismatched sizes");
//    if(A.size() > std::numeric_limits<LAPACK_INT>::max()) 
//        throw std::runtime_error("vecref*vecref: overflow of size beyond LAPACK_INT range");
//#endif
//    return ddot_wrapper(A.size(),A.cstore(),A.stride(),B.cstore(),B.stride());
//    }
//
//vec
//randomVec(long size)
//    {
//    vec res(size);
//    res.randomize();
//    return res;
//    }

//bool
//overlaps(const Real* b1, const Real* e1,
//         const Real* b2, const Real* e2)
//    {
//    using ptr_less = std::less<const Real*>;
//    return (ptr_less(b1,b2) && ptr_less(b2,e1))
//        || (ptr_less(b1,e2) && ptr_less(e2,e1));
//    }


}; //namespace itensor
