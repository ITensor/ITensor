//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SIMPLEMATRIX_H_
#define __ITENSOR_SIMPLEMATRIX_H_

#include "types.h"
#include "print.h"

namespace itensor {

class vec;

class vecref
    {
    public:
    using iterator = Real*;
    using const_iterator = const Real*;
    using value_type = Real;
    using size_type = long;
    private:
    Real *store_ = nullptr;
    const Real *cstore_ = nullptr;
    size_type size_ = 0;
    public:

    vecref() { }

    vecref(long size) : size_(size) { }

    vecref(const Real* sto, 
           long size)
        :
        store_(nullptr),
        cstore_(sto),
        size_(size)
        { }

    vecref(Real* sto, 
           long size)
        :
        store_(sto),
        cstore_(sto),
        size_(size)
        { }

    void virtual
    operator=(const vecref& other);

    size_type
    size() const { return size_; }

    explicit operator bool() const { return bool(cstore_); }

    bool
    readOnly() const { return !bool(store_); }

    const Real*
    store() const { return cstore_; }
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
    operator()(long i) const { return cstore_[i-1]; }
    Real&
    operator()(long i)
        { 
#ifdef DEBUG
        if(readOnly()) throw std::runtime_error("vecref is read only");
#endif
        return store_[i-1];
        }

    iterator
    begin() 
        { 
#ifdef DEBUG
        if(readOnly()) throw std::runtime_error("vecref is read only");
#endif
        return store_; 
        }
    iterator
    end() 
        { 
#ifdef DEBUG
        if(readOnly()) throw std::runtime_error("vecref is read only");
#endif
        return store_+size_; 
        }
    const_iterator
    begin() const{ return cstore_; }
    const_iterator
    end() const{ return cstore_+size_; }
    const_iterator
    cbegin() const{ return cstore_; }
    const_iterator
    cend() const{ return cstore_+size_; }

    vecref(const vec& other) = delete;
    vecref&
    operator=(const vec& other) = delete;
    };

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

    private:

    void
    assignFromRef(const vecref& other)
        {
        parent::operator=(other);
        data_ = storage_type(other.cbegin(),other.cend());
        store(data_.data());
        }

    void
    assignFromVec(const vec& other)
        {
        const vecref& oref = other;
        parent::operator=(oref);
        data_ = other.data_;
        }

    void
    moveFromVec(vec&& other)
        {
        const vecref& oref = other;
        parent::operator=(oref);
        data_ = std::move(other.data_);
        }
    public:
    const Real*
    data() const { return data_.data(); }
    };

class matrix;

class matrixref
    {
    public:
    using iterator = Real*;
    using const_iterator = const Real*;
    using value_type = Real;
    using size_type = long;
    private:
    Real *store_ = nullptr;
    const Real *cstore_ = nullptr;
    long nrows_ = 0, 
         ncols_ = 0;
    bool transpose_ = false;
    public:

    matrixref() { }

    matrixref(long nro, 
              long ncol, 
              bool trans = false);

    matrixref(const Real* sto, 
              long nro, 
              long ncol, 
              bool trans = false);

    matrixref(Real* sto, 
              long nro, 
              long ncol, 
              bool trans = false);

    void virtual
    operator=(const matrixref& other);

    matrixref(const matrix& other) = delete;
    matrixref&
    operator=(const matrix& other) = delete;

    long
    Nrows() const { return transpose_ ? ncols_ : nrows_; }
    long
    Ncols() const { return transpose_ ? nrows_ : ncols_; }

    explicit operator bool() const { return bool(cstore_); }

    bool
    transpose() const { return transpose_; }
    void
    applyTrans() { transpose_ = !transpose_; }

    bool
    readOnly() const { return !bool(store_); }

    const Real*
    store() const { return cstore_; }
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

    matrixref 
    t();

    Real
    operator()(long i, long j) const { return cstore_[index(i,j)]; }
    Real&
    operator()(long i, long j) 
        { 
#ifdef DEBUG
        if(readOnly()) throw std::runtime_error("matrixref is read only");
#endif
        return store_[index(i,j)]; 
        }

    long
    size() const { return long(nrows_*ncols_); }

    iterator
    begin() 
        { 
#ifdef DEBUG
        if(readOnly()) throw std::runtime_error("matrixref is read only");
#endif
        return store_; 
        }
    iterator
    end() 
        { 
#ifdef DEBUG
        if(readOnly()) throw std::runtime_error("matrixref is read only");
#endif
        return store_+size(); 
        }
    const_iterator
    begin() const{ return cstore_; }
    const_iterator
    end() const{ return cstore_+size(); }
    const_iterator
    cbegin() const{ return cstore_; }
    const_iterator
    cend() const{ return cstore_+size(); }

    private:
    
    long 
    index0(long i, long j) const
        { 
        return transpose_ ? j+i*nrows_ : i+j*nrows_;
        }
    long 
    index(long i, long j) const { return index0(i-1,j-1); }
    };

class matrix : public matrixref
    {
    public:
    using parent = matrixref;
    using storage_type = std::vector<Real>;
    using iterator = parent::iterator;
    using const_iterator = parent::const_iterator;
    using value_type = parent::value_type;
    using size_type = parent::size_type;
    private:
    storage_type data_;
    public:

    matrix() { }

    matrix(long nro, 
           long ncol, 
           bool trans = false)
        : matrixref(nro,ncol,trans)
        {
        data_ = storage_type(nro*ncol,0);
        store(data_.data());
        }
      
    matrix(const matrix& other) { assignFrom(other); }

    matrix(matrix&& other) { moveFrom(std::move(other)); }

    matrix(const matrixref& other) { assignFromRef(other); }

    matrix& 
    operator=(const matrix& other) { assignFrom(other); return *this; }
    matrix& 
    operator=(matrix&& other) { moveFrom(std::move(other)); return *this; }
    void virtual
    operator=(const matrixref& other) override { assignFromRef(other); }

    private:
    void
    assignFromRef(const matrixref& other)
        {
        parent::operator=(other);
        data_ = storage_type(other.cbegin(),other.cend());
        store(data_.data());
        }

    void
    assignFrom(const matrix& other)
        {
        const matrixref& oref = other;
        parent::operator=(oref);
        data_ = other.data_;
        }

    void
    moveFrom(matrix&& other)
        {
        const matrixref& oref = other;
        parent::operator=(oref);
        data_ = std::move(other.data_);
        }
    public:
    const Real*
    data() const { return data_.data(); }
    };

std::ostream&
operator<<(std::ostream& s, const vecref& v);
std::ostream&
operator<<(std::ostream& s, const matrixref& M);

// C += A*B
void
mult_add(const matrixref& A, 
         const matrixref& B, 
         matrixref& C);

// C = A*B
void
mult(const matrixref& A, 
     const matrixref& B, 
     matrixref& C);

};

#endif
