//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MATRIX___H_
#define __ITENSOR_MATRIX___H_

#include "vector.h"

namespace itensor {

class matrix;

struct mrange
    {
    long rn=0,rs=0,cn=0,cs=0;
    mrange() { } 
    mrange(long rn_, long rs_,
         long cn_, long cs_) 
        : rn(rn_),rs(rs_),cn(cn_),cs(cs_) 
        { }
    mrange(long rn_, long cn_)
        : rn(rn_),rs(1),cn(cn_),cs(rn_) 
        { }
    long
    index(long i, long j) const { return (i-1)*rs+(j-1)*cs; }
    long
    index0(long i, long j) const { return i*rs+j*cs; }
    long
    area() const { return rn*cn; }
    };

mrange inline
transpose(const mrange& ind)
    {
    return mrange(ind.cn,ind.cs,ind.rn,ind.rs); 
    }

class matrixref
    {
    public:
    using iterator = Real*;
    using const_iterator = const Real*;
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

    explicit operator bool() const { return bool(cstore_); }

    void
    applyTrans() { ind_ = transpose(ind_); }

    bool
    readOnly() const { return !bool(store_); }

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

    long
    size() const { return ind_.area(); }

//    iterator
//    begin() 
//        { 
//#ifdef DEBUG
//        if(readOnly()) throw std::runtime_error("matrixref is read only");
//#endif
//        return store_; 
//        }
//    iterator
//    end() 
//        { 
//#ifdef DEBUG
//        if(readOnly()) throw std::runtime_error("matrixref is read only");
//#endif
//        return store_+size(); 
//        }
//    const_iterator
//    begin() const { return cstore_; }
//    const_iterator
//    end() const { return cstore_+size(); }
//    const_iterator
//    cbegin() const { return cstore_; }
//    const_iterator
//    cend() const { return cstore_+size(); }

    };

vecref
diagonal(const matrixref& m) 
    { 
    auto vsize = std::min(m.Nrows(),m.Ncols());
    auto vstrd = m.rowStride()+m.colStride();
    if(m.readOnly()) return vecref(m.cstore(),vsize,vstrd);
    else           return vecref(m.store(),vsize,vstrd);
    }

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
        throw std::runtime_error("assignFromRef currently broken");
        //data_ = storage_type(other.cbegin(),other.cend());

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
