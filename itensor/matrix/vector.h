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
           long stride = 1)
        :
        store_(nullptr),
        cstore_(sto),
        strd_(stride),
        size_(size)
        { }

    vecref(Real* sto, 
           long size,
           long stride = 1)
        :
        store_(sto),
        cstore_(sto),
        strd_(stride),
        size_(size)
        { }

    void virtual
    operator=(const vecref& other);

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
    store() const
        { 
#ifdef DEBUG
        if(readOnly()) throw std::runtime_error("vecref read-only: call cstore() or call store() on const vecref object");
#endif
        return store_; 
        }

    void
    store(const Real* newstore) 
        { 
        store_ = nullptr;
        cstore_ = newstore;
        }
    void
    store(Real* newstore) 
        { 
        store_ = newstore;
        cstore_ = newstore;
        }

    Real
    operator()(long i) const { return cstore_[(i-1)*strd_]; }
    Real&
    operator()(long i)
        { 
#ifdef DEBUG
        if(readOnly()) throw std::runtime_error("vecref is read only");
#endif
        return store_[(i-1)*strd_];
        }
    Real
    get(long i) const { return cstore_[(i-1)*strd_]; }

    iterator
    begin() 
        { 
#ifdef DEBUG
        if(readOnly()) throw std::runtime_error("vecref is read only");
#endif
        return iterator(store_,strd_); 
        }
    iterator
    end() 
        { 
#ifdef DEBUG
        if(readOnly()) throw std::runtime_error("vecref is read only");
#endif
        return iterator(store_+size_*strd_,strd_); 
        }
    const_iterator
    begin() const{ return const_iterator(cstore_,strd_); }
    const_iterator
    end() const{ return const_iterator(cstore_+size_*strd_,strd_); }
    const_iterator
    cbegin() const{ return const_iterator(cstore_,strd_); }
    const_iterator
    cend() const{ return const_iterator(cstore_+size_*strd_,strd_); }

    void virtual
    clear() { *this = vecref(); }

    vecref(const vec& other) = delete;
    vecref&
    operator=(const vec& other) = delete;


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
    assignFromRef(const vecref& other)
        {
        if(&other == this) return;
        data_ = storage_type(other.cbegin(),other.cend());
        store(data_.data());
        parent::size(data_.size());
        parent::stride(1);
        }

    void
    assignFromVec(const vec& other)
        {
        if(&other == this) return;
        data_ = other.data_;
        store(data_.data());
        parent::size(data_.size());
        parent::stride(1);
        }

    void
    moveFromVec(vec&& other)
        {
        const vecref& oref = other;
        parent::operator=(oref);
        data_ = std::move(other.data_);
        other.clear();
        }
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
