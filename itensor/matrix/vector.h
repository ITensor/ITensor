//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_VECTOR___H_
#define __ITENSOR_VECTOR___H_

#include "types.h"
#include "print.h"
#include "strideiter.h"

namespace itensor {

class vec;

class vecref
    {
    public:
    using iterator = stride_iter<Real*>;
    using const_iterator = stride_iter<const Real*>;
    using value_type = Real;
    using size_type = long;
    private:
    Real *store_ = nullptr;
    const Real *cstore_ = nullptr;
    size_type strd_ = 1;
    size_type size_ = 0;
    public:

    vecref() { }

    vecref(long size) : size_(size) { }

    vecref(const Real* sto, 
           long size,
           long stride = 1);

    vecref(Real* sto, 
           long size,
           long stride = 1);

    void virtual
    operator=(const vecref& other);

    vecref(const vec& other) = delete;

    vecref&
    operator=(const vec& other) { assignFromVec(other); return *this; }

    size_type
    size() const { return size_; }
    size_type
    stride() const { return strd_; }
    bool
    contiguous() const { return strd_ == 1; }
    bool
    readOnly() const { return !bool(store_); }

    explicit operator bool() const { return bool(cstore_); }

    void
    operator*=(Real fac);
    void
    operator/=(Real fac);
    void
    operator+=(const vecref& other);
    void
    operator-=(const vecref& other);

    const Real*
    cstore() const { return cstore_; }

    Real*
    store() const;

    void
    store(const Real* newstore);
    void
    store(Real* newstore);

    Real
    operator()(long i) const;
    Real&
    operator()(long i);
    Real
    get(long i) const { return cstore_[(i-1)*strd_]; }

    iterator
    begin();
    iterator
    end();
    const_iterator
    begin() const;
    const_iterator
    end() const;
    const_iterator
    cbegin() const;
    const_iterator
    cend() const;

    void virtual
    clear() { *this = vecref(); }

    private:

    void virtual
    assignFromVec(const vec& other);

    public:
    void
    size(size_type nsize) { size_ = nsize; }
    void
    stride(size_type nstride) { strd_ = nstride; }
    };

vecref inline
subVector(const vecref& v,
          long start,
          long stop)
    {
    if(v.readOnly()) return vecref(v.cstore()+(start-1),stop-start+1);
    else             return vecref(v.store()+(start-1),stop-start+1);
    }

class vec : public vecref
    {
    public:
    using parent = vecref;
    using storage_type = std::vector<Real>;
    using iterator = parent::iterator;
    using const_iterator = parent::const_iterator;
    using value_type = parent::value_type;
    using size_type = parent::size_type;
    public:
    storage_type data_;
    public:

    vec() { }

    vec(long size) : parent(size)
        {
        data_ = storage_type(size);
        store(data_.data());
        }

    vec(const vec& other) { assignFromVec(other); }

    vec(vec&& other) { moveFromVec(std::move(other)); }

    vec(const vecref& other) { assignFromRef(other); }

    vec&
    operator=(const vec& other) { assignFromVec(other); return *this; }
    vec& 
    operator=(vec&& other) { moveFromVec(std::move(other)); return *this; }
    void virtual
    operator=(const vecref& other) override { assignFromRef(other); }

    void virtual
    clear() override
        {
        parent::clear();
        data_.clear();
        }

    private:

    void
    assignFromRef(const vecref& other);

    void virtual
    assignFromVec(const vec& other) override;

    void
    moveFromVec(vec&& other);

    public:
    const Real*
    data() const { return data_.data(); }
    };

vec inline
operator+(vec A, const vec& B)
    {
    A += B;
    return A;
    }
vec inline
operator+(const vec& A, vec&& B)
    {
    vec res(std::move(B));
    res += A;
    return res;
    }
vec inline
operator-(vec A, const vec& B)
    {
    A -= B;
    return A;
    }
vec inline
operator-(const vec& A, vec&& B)
    {
    vec res(std::move(B));
    res *= -1;
    res += A;
    return res;
    }

Real
norm(const vecref& v);

vec
randomVec(long size);

std::ostream&
operator<<(std::ostream& s, const vecref& v);

};

#endif
