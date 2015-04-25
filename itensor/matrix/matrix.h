//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MATRIX___H_
#define __ITENSOR_MATRIX___H_

#include "vector.h"
#include "miterator.h"

namespace itensor {

class matrix;

class matrixref
    {
    public:
    using iterator = miterator<Real*>;
    using const_iterator = miterator<const Real*>;
    using value_type = Real;
    using size_type = long;
    private:
    mrange ind_;
    Real *store_ = nullptr;
    const Real *cstore_ = nullptr;
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

    matrixref(Real* sto, 
              const mrange& ind);
    matrixref(const Real* sto, 
              const mrange& ind);

    void virtual
    operator=(const matrixref& other);

    matrixref(const matrix& other) = delete;

    matrixref&
    operator=(const matrix& other) = delete;

    long
    Nrows() const { return ind_.rn; }
    long
    Ncols() const { return ind_.cn; }
    long
    rowStride() const { return ind_.rs; }
    long
    colStride() const { return ind_.cs; }

    bool
    transposed() const { return isTransposed(ind_); }
    bool
    contiguous() const { return isContiguous(ind_); }
    bool
    readOnly() const { return !bool(store_); }

    explicit operator bool() const { return bool(cstore_); }

    void
    applyTrans() { ind_ = transpose(ind_); }
    matrixref 
    t();
    void
    randomize();

    void
    operator*=(Real fac);
    void
    operator/=(Real fac);

    void
    operator+=(const matrixref& other);
    void
    operator-=(const matrixref& other);

    Real
    operator()(long i, long j) const { return cstore_[ind_.index(i,j)]; }
    Real&
    operator()(long i, long j) 
        { 
#ifdef DEBUG
        if(readOnly()) throw std::runtime_error("matrixref is read only");
#endif
        return store_[ind_.index(i,j)]; 
        }
    Real
    get(long i, long j) const { return cstore_[ind_.index(i,j)];  }

    long
    size() const { return ind_.area(); }

    const mrange&
    ind() const { return ind_; }
    void
    ind(const mrange& ni) { ind_ = ni; }

    const Real*
    cstore() const { return cstore_; }

    Real*
    store() const
        { 
#ifdef DEBUG
        if(readOnly()) throw std::runtime_error("matrixref read-only: call cstore() or call store() on const matrixref object");
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

    iterator
    begin() 
        { 
#ifdef DEBUG
        if(readOnly()) throw std::runtime_error("matrixref is read only");
#endif
        return iterator(store_,ind_); 
        }
    iterator
    end() 
        { 
#ifdef DEBUG
        if(readOnly()) throw std::runtime_error("matrixref is read only");
#endif
        return iterator(ind_);
        }
    const_iterator
    begin() const { return const_iterator(cstore_,ind_); }
    const_iterator
    end() const { return const_iterator(ind_); }
    const_iterator
    cbegin() const { return const_iterator(cstore_,ind_); }
    const_iterator
    cend() const { return const_iterator(ind_); }
    void virtual
    clear() { *this = matrixref(); }
    };

vecref
diagonal(const matrixref& m);

vecref
column(const matrixref& m, long j);

vecref
row(const matrixref& m, long j);

matrixref
subMatrix(const matrixref& m,
          long rstart,
          long rstop,
          long cstart,
          long cstop);


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

// y = M*x
void
mult(const matrixref& M,
     const vecref& x, 
     vec& y,
     bool fromleft = false);

void
diagSymmetric(const matrixref& M,
              matrixref& U,
              vecref& d);

void
diagSymmetric(const matrixref& M,
              matrix& U,
              vec& d);

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

    void
    operator+=(const matrix& other);
    void
    operator-=(const matrix& other);

    void virtual
    clear() override
        {
        parent::clear();
        data_.clear();
        }

    private:
    void
    assignFromRef(const matrixref& other)
        {
        if(&other == this) return;
        data_ = storage_type(other.cbegin(),other.cend());
        store(data_.data());
        parent::ind(mrange(other.Nrows(),other.Ncols()));
        }

    void
    assignFrom(const matrix& other)
        {
        if(&other == this) return;
        data_ = other.data_;
        store(data_.data());
        parent::ind(mrange(other.Nrows(),other.Ncols()));
        }

    void
    moveFrom(matrix&& other)
        {
        const matrixref& oref = other;
        parent::operator=(oref);
        data_ = std::move(other.data_);
        other.clear();
        }

    public:
    const Real*
    data() const { return data_.data(); }
    };

matrix inline
operator*(const matrixref& A,
          const matrixref& B)
    {
    matrix C(A.Nrows(),B.Ncols());
    mult(A,B,C);
    return C;
    }

vec inline
operator*(const matrixref& M,
          const vecref& v)
    {
    vec res(M.Nrows());
    mult(M,v,res);
    return res;
    }

vec inline
operator*(const vecref& v,
          const matrixref& M)
    {
    vec res(M.Ncols());
    mult(M,v,res,true);
    return res;
    }

matrix inline
operator+(matrix A, const matrix& B)
    {
    A += B;
    return A;
    }
matrix inline
operator+(const matrix& A, matrix&& B)
    {
    matrix res(std::move(B));
    res += A;
    return res;
    }
matrix inline
operator-(matrix A, const matrix& B)
    {
    A -= B;
    return A;
    }
matrix inline
operator-(const matrix& A, matrix&& B)
    {
    matrix res(std::move(B));
    res *= -1;
    res += A;
    return res;
    }
matrix inline
operator*(matrix A, Real fac) { A *= fac; return A; }
matrix inline
operator*(Real fac, matrix A) { return operator*(A,fac); }
matrix inline
operator/(matrix A, Real fac) { A /= fac; return A; }

matrix
randomMatrix(long Nr, long Nc);

Real
norm(const matrixref& M);

Real
norm(const matrix& M);

std::ostream&
operator<<(std::ostream& s, const matrixref& M);

};

#endif
