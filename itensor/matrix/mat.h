//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MAT__H_
#define __ITENSOR_MAT__H_

#include "vec.h"
#include "miterator.h"

namespace itensor {

template<typename T>
class MatRefT;

using MatRef = MatRefT<Real>;

using CMatRef = MatRefT<const Real>;

template<typename T>
class MatRefT
    {
    public:
    using iterator = miterator<T*>;
    using const_iterator = miterator<const T*>;
    using value_type = std::remove_const_t<T>;
    using pointer = T*;
    using reference = T&;
    using size_type = long;
    private:
    pointer pdata_ = nullptr;
    mrange ind_;
    public:

    MatRefT() { }

    MatRefT(pointer pdata, 
            long nrows,
            long ncols,
            bool trans = false)
        :
        pdata_(pdata),
        ind_(trans ? mrange(ncols,nrows,nrows,1) : mrange(nrows,ncols))
        { }

    MatRefT(pointer pdata, 
            const mrange& ind)
        :
        pdata_(pdata),
        ind_(ind)
        { }

    MatRefT(const MatRefT<value_type>& other)
        :
        pdata_(other.pdata_),
        ind_(other.ind_)
        { }

    MatRefT(const MatRefT<const value_type>& other)
        :
        pdata_(other.pdata_),
        ind_(other.ind_)
        { }

    MatRefT&
    operator=(const MatRefT<value_type>& other)
        {
        pdata_ = other.pdata_;
        ind_ = other.ind_;
        return *this;
        }

    MatRefT&
    operator=(const MatRefT<const value_type>& other)
        {
        pdata_ = other.pdata_;
        ind_ = other.ind_;
        return *this;
        }

    long
    Nrows() const { return ind_.rn; }
    long
    Ncols() const { return ind_.cn; }
    long
    rowStride() const { return ind_.rs; }
    long
    colStride() const { return ind_.rs; }

    size_type
    size() const { return ind_.area(); }

    const mrange&
    ind() const { return ind_; }

    bool
    contiguous() const { return isContiguous(ind_); }

    bool
    transposed() const { return isTransposed(ind_); }

    explicit operator bool() const { return bool(pdata_); }

    pointer
    data() const { return pdata_; }

    void
    applyTrans() { ind_ = transpose(ind_); }

    MatRefT 
    t() const
        {
        MatRefT res(*this);
        res.applyTrans();
        return res;
        }

    reference
    operator()(long i, long j) const { return pdata_[ind_.index(i,j)]; }

    iterator 
    begin() const { return iterator(pdata_,ind_); }

    iterator 
    end() const { return iterator(ind_); }

    const_iterator 
    cbegin() const { return const_iterator(pdata_,ind_); }

    const_iterator 
    cend() const { return const_iterator(ind_); }

    friend MatRefT<value_type>;
    friend MatRefT<const value_type>;
    };

std::ostream&
operator<<(std::ostream& s, CMatRef M);

MatRef inline
makeMatRef(Real* p, 
           long nrows,
           long ncols,
           bool trans = false)
    {
    return MatRef(p,nrows,ncols,trans);
    }

CMatRef inline
makeMatRef(const Real* cp, 
           long nrows,
           long ncols,
           bool trans = false)
    {
    return CMatRef(cp,nrows,ncols,trans); 
    }


//Copy data referenced by b to memory referenced by a
MatRef&
operator&=(MatRef& a, CMatRef b);

MatRef&
operator*=(MatRef& v, Real fac);

MatRef&
operator/=(MatRef& v, Real fac);

MatRef&
operator+=(MatRef& a, CMatRef b);

MatRef&
operator-=(MatRef& a, CMatRef b);

// compute matrix multiply (dgemm) A*B
// write result to memory referenced by C
void
mult(CMatRef A, 
     CMatRef B, 
     MatRef  C);

};

#endif
