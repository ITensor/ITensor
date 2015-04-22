//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SIMPLEMATRIX_H_
#define __ITENSOR_SIMPLEMATRIX_H_

#include "types.h"
#include "print.h"

namespace itensor {

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

    matrixref(const matrixref& other) = default;

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
    std::vector<Real> data_;
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
      
    matrix(const matrix& other) = default;
    matrix(matrix&& other) = default;

    matrix(const matrixref& other)
        :
        matrixref(other)
        { 
        data_ = storage_type(other.cbegin(),other.cend());
        store(data_.data());
        }

    };

std::ostream&
operator<<(std::ostream& s, const matrixref& M);

// C += A*B
void
mult_add(matrixref A, 
         matrixref B, 
         matrixref C);

// C = A*B
void
mult(matrixref A, 
     matrixref B, 
     matrixref C);

};

#endif
