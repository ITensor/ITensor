//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_VEC_H_
#define __ITENSOR_VEC_H_

#include "types.h"
#include "print.h"
#include "strideiter.h"

namespace itensor {

template<typename T>
class VecRefT;

using VecRef = VecRefT<Real>;

using CVecRef = VecRefT<const Real>;

template<typename T>
class VecRefT
    {
    public:
    using iterator = stride_iter<T*>;
    using const_iterator = stride_iter<const T*>;
    using value_type = std::remove_const_t<T>;
    using pointer = T*;
    using reference = T&;
    using size_type = long;
    private:
    pointer pdata_ = nullptr;
    size_type strd_ = 1;
    size_type size_ = 0;
    public:

    VecRefT() { }

    VecRefT(pointer pdata, 
            long size,
            long stride = 1)
        :
        pdata_(pdata),
        strd_(stride),
        size_(size)
        { }

    operator VecRefT<const T>() const { return VecRefT<const T>(pdata_,size_,strd_); }

    size_type
    size() const { return size_; }

    size_type
    stride() const { return strd_; }

    bool
    contiguous() const { return strd_ == 1; }

    explicit operator bool() const { return bool(pdata_); }

    pointer
    data() const { return pdata_; }

    reference
    operator()(long i) const { return pdata_[(i-1)*strd_]; }

    iterator 
    begin() const { return iterator(pdata_,strd_); }

    iterator 
    end() const { return iterator(pdata_+size_*strd_,strd_); }

    const_iterator 
    cbegin() const { return const_iterator(pdata_,strd_); }

    const_iterator 
    cend() const { return const_iterator(pdata_+size_*strd_,strd_); }
    };

VecRef inline
makeVecRef(Real* pd, 
        long size,
        long stride = 1)
    {
    return VecRef(pd,size,stride); 
    }

CVecRef inline
makeVecRef(const Real* cpd, 
        long size,
        long stride = 1)
    {
    return CVecRef(cpd,size,stride); 
    }

//Copy data referenced by b to memory referenced by a
VecRef&
operator&=(VecRef& a, CVecRef b);

VecRef&
operator*=(VecRef& v, Real fac);

VecRef&
operator/=(VecRef& v, Real fac);

VecRef&
operator+=(VecRef& a, CVecRef b);

VecRef&
operator-=(VecRef& a, CVecRef b);

//Dot product
Real
operator*(CVecRef a, CVecRef b);

class Vec
    {
    public:
    using storage_type = std::vector<Real>;
    using iterator = storage_type::iterator;
    using const_iterator = storage_type::const_iterator;
    using value_type = Real;
    using size_type = storage_type::size_type;
    public:
    storage_type data_;
    public:

    Vec() { }

    explicit
    Vec(long size) : data_(size) { }

    Vec(const Vec& other) { assignFromVec(other); }

    Vec(Vec&& other) { moveFromVec(std::move(other)); }

    explicit
    Vec(CVecRef ref) { assignFromRef(ref); }

    Vec&
    operator=(const Vec& other) { assignFromVec(other); return *this; }
    Vec& 
    operator=(Vec&& other) { moveFromVec(std::move(other)); return *this; }
    Vec&
    operator=(CVecRef ref) { assignFromRef(ref); return *this; }

    operator CVecRef() const { return CVecRef(data_.data(),data_.size()); }

    Real&
    operator()(long i) 
        { 
#ifdef DEBUG
        return data_.at(i-1); 
#else
        return data_[i-1]; 
#endif
        }

    Real
    operator()(long i) const 
        { 
#ifdef DEBUG
        return data_.at(i-1); 
#else
        return data_[i-1]; 
#endif
        }

    Vec&
    operator*=(Real fac);
    Vec&
    operator/=(Real fac);
    Vec&
    operator+=(const Vec& other);
    Vec&
    operator-=(const Vec& other);
    Vec&
    operator+=(CVecRef other);
    Vec&
    operator-=(CVecRef other);

    explicit operator bool() const { return !data_.empty(); }

    Real*
    data() { return data_.data(); }

    const Real*
    data() const { return data_.data(); }

    size_type
    size() const { return data_.size(); }

    void
    clear() { data_.clear(); }

    iterator
    begin() { return data_.begin(); }
    iterator
    end() { return data_.end(); }
    const_iterator
    begin() const { return data_.cbegin(); }
    const_iterator
    end() const { return data_.cend(); }
    const_iterator
    cbegin() const { return data_.cbegin(); }
    const_iterator
    cend() const { return data_.cend(); }

    private:

    void
    assignFromRef(CVecRef other)
        {
        //Copy data from other contiguously into data_
        data_ = storage_type(other.cbegin(),other.cend());
        }

    void
    assignFromVec(const Vec& other)
        {
        if(&other == this) return;
        data_ = other.data_;
        }

    void
    moveFromVec(Vec&& other)
        {
        data_ = std::move(other.data_);
        other.clear();
        }

    };

VecRef inline
makeVecRef(Vec& v, 
        long size,
        long stride = 1)
    {
    return VecRef(v.data(),size,stride); 
    }

VecRef inline
makeVecRef(Vec& v)
    {
    return VecRef(v.data(),v.size()); 
    }

CVecRef inline
makeVecRef(const Vec& v, 
        long size,
        long stride = 1)
    {
    return CVecRef(v.data(),size,stride); 
    }

CVecRef inline
makeVecRef(const Vec& v)
    {
    return CVecRef(v.data(),v.size()); 
    }

Vec inline
operator+(Vec A, const Vec& B) { A += B; return A; }
Vec inline
operator+(const Vec& A, Vec&& B) { Vec res(std::move(B)); res += A; return res; }
Vec inline
operator-(Vec A, const Vec& B) { A -= B; return A; }
Vec inline
operator-(const Vec& A, Vec&& B) { Vec res(std::move(B)); res *= -1; res += A; return res; }
inline Vec& Vec::
operator*=(Real fac) { auto r = makeVecRef(*this); r *= fac; return *this; }
inline Vec& Vec::
operator/=(Real fac) { auto r = makeVecRef(*this); r /= fac; return *this; }
inline Vec& Vec::
operator+=(const Vec& other) { auto r = makeVecRef(*this); r += makeVecRef(other); return *this; }
inline Vec& Vec::
operator-=(const Vec& other) { auto r = makeVecRef(*this); r -= makeVecRef(other); return *this; }
inline Vec& Vec::
operator+=(CVecRef other) { auto r = makeVecRef(*this); r += other; return *this; }
inline Vec& Vec::
operator-=(CVecRef other) { auto r = makeVecRef(*this); r -= other; return *this; }

//Copy contents of Vec to memory referenced by VecRef
inline VecRef&
operator&=(VecRef& ref, const Vec& v) { return operator&=(ref,makeVecRef(v)); }

//Dot product
Real inline
operator*(const Vec& a, const Vec& b) { return operator*(makeVecRef(a),makeVecRef(b)); }

Real
norm(CVecRef v);

Real inline
norm(const Vec& v) { return norm(makeVecRef(v)); }

VecRef
randomize(VecRef v);

inline Vec&
randomize(Vec& v) { randomize(makeVecRef(v)); return v; }

Vec
randomVec(long size);

std::ostream&
operator<<(std::ostream& s, CVecRef v);

inline std::ostream&
operator<<(std::ostream& s, const Vec& v) { return operator<<(s,makeVecRef(v)); }

//Return ref to elements [start,stop] of a Vec, inclusive
template<typename Vec_>
auto
subVector(Vec_& v,
          long start,
          long stop)
    {
    return makeVecRef(v.data()+(start-1),stop-start+1);
    }

};

#endif
