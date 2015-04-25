//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "vector.h"
#include "lapack_wrap.h"
#include <limits>
#include "detail/algs.h"

namespace itensor {

vecref::
vecref(const Real* sto, 
       long size,
       long stride)
    :
    store_(nullptr),
    cstore_(sto),
    strd_(stride),
    size_(size)
    { }

vecref::
vecref(Real* sto, 
       long size,
       long stride)
    :
    store_(sto),
    cstore_(sto),
    strd_(stride),
    size_(size)
    { }

void vecref::
operator=(const vecref& other)
    {
    store_ = other.store_;
    cstore_ = other.cstore_;
    strd_ = other.strd_;
    size_ = other.size_;
    }

void vecref::
assignFromVec(const vec& other)
    {
    if(&other == this) return;
    if(readOnly()) throw std::runtime_error("vecref is read only, cannot assign from vec");
#ifdef DEBUG
    if(other.size() != size_) throw std::runtime_error("mismatched sizes in vecref::operator=(vec)");
#endif
    auto po = other.cbegin();
    for(auto& el : *this) 
        {
        el = *po;
        ++po;
        }
    }

Real* vecref::
store() const
    { 
#ifdef DEBUG
    if(readOnly()) throw std::runtime_error("vecref read-only: call cstore() or call store() on const vecref object");
#endif
    return store_; 
    }

void vecref::
store(const Real* newstore) 
    { 
    store_ = nullptr;
    cstore_ = newstore;
    }

void vecref::
store(Real* newstore) 
    { 
    store_ = newstore;
    cstore_ = newstore;
    }

Real vecref::
operator()(long i) const { return cstore_[(i-1)*strd_]; }

Real& vecref::
operator()(long i)
    { 
#ifdef DEBUG
    if(readOnly()) throw std::runtime_error("vecref is read only");
#endif
    return store_[(i-1)*strd_];
    }

void vecref::
randomize()
    {
#ifdef DEBUG
    if(readOnly()) throw std::runtime_error("randomize: vecref is read only");
#endif
    for(auto& el : *this) el = detail::quickran();
    }

vecref::iterator vecref::
begin() 
    { 
#ifdef DEBUG
    if(readOnly()) throw std::runtime_error("vecref is read only");
#endif
    return iterator(store_,strd_); 
    }
vecref::iterator vecref::
end() 
    { 
#ifdef DEBUG
    if(readOnly()) throw std::runtime_error("vecref is read only");
#endif
    return iterator(store_+size_*strd_,strd_); 
    }

vecref::const_iterator vecref::
begin() const { return const_iterator(cstore_,strd_); }

vecref::const_iterator vecref::
end() const { return const_iterator(cstore_+size_*strd_,strd_); }

vecref::const_iterator vecref::
cbegin() const { return const_iterator(cstore_,strd_); }

vecref::const_iterator vecref::
cend() const { return const_iterator(cstore_+size_*strd_,strd_); }


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

void vec::
assignFromRef(const vecref& other)
    {
    if(&other == this) return;
    data_ = storage_type(other.cbegin(),other.cend());
    store(data_.data());
    parent::size(data_.size());
    parent::stride(1);
    }

void vec::
assignFromVec(const vec& other)
    {
    if(&other == this) return;
    data_ = other.data_;
    store(data_.data());
    parent::size(data_.size());
    parent::stride(1);
    }

void vec::
moveFromVec(vec&& other)
    {
    const vecref& oref = other;
    parent::operator=(oref);
    data_ = std::move(other.data_);
    other.clear();
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

Real
operator*(const vecref& A, const vecref& B)
    {
#ifdef DEBUG
    if(A.size() != B.size()) throw std::runtime_error("vecref*vecref: mismatched sizes");
    if(A.size() > std::numeric_limits<LAPACK_INT>::max()) 
        throw std::runtime_error("vecref*vecref: overflow of size beyond LAPACK_INT range");
#endif
    return ddot_wrapper(A.size(),A.cstore(),A.stride(),B.cstore(),B.stride());
    }

vec
randomVec(long size)
    {
    vec res(size);
    res.randomize();
    return res;
    }

}; //namespace itensor
