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
    transposed() const { return (ind_.rs==ind_.cn && ind_.cs==1); }
    bool
    contiguous() const { return (ind_.rs==1 && ind_.cs==ind_.rn) || transposed(); }

    explicit operator bool() const { return bool(cstore_); }

    void
    applyTrans() { ind_ = transpose(ind_); }

    bool
    readOnly() const { return !bool(store_); }


    matrixref 
    t();

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
        return make_end(iterator(store_,ind_));
        }
    const_iterator
    begin() const { return const_iterator(cstore_,ind_); }
    const_iterator
    end() const { return make_end(const_iterator(cstore_,ind_)); }
    const_iterator
    cbegin() const { return const_iterator(cstore_,ind_); }
    const_iterator
    cend() const { return make_end(const_iterator(cstore_,ind_)); }
    };

vecref
diagonal(const matrixref& m);

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

void
diagHermitian(const matrixref& M,
              matrixref& U,
              vecref& d);

void
diagHermitian(const matrixref& M,
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

matrix inline
operator*(const matrixref& A,
          const matrixref& B)
    {
    matrix C(A.Nrows(),B.Ncols());
    mult(A,B,C);
    return C;
    }


std::ostream&
operator<<(std::ostream& s, const vecref& v);
std::ostream&
operator<<(std::ostream& s, const matrixref& M);


};

#endif
