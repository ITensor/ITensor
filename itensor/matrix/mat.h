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
    MRange ind_;
    public:

    MatRefT() { }

    MatRefT(pointer pdata, 
            long nrows,
            long ncols,
            bool trans = false)
        :
        pdata_(pdata),
        ind_(trans ? MRange(ncols,nrows,nrows,1) : MRange(nrows,ncols))
        { }

    MatRefT(pointer pdata, 
            const MRange& ind)
        :
        pdata_(pdata),
        ind_(ind)
        { }

    operator MatRefT<const T>() const { return MatRefT<const T>(pdata_,ind_); }

    long
    Nrows() const { return ind_.rn; }
    long
    Ncols() const { return ind_.cn; }
    long
    rowStride() const { return ind_.rs; }
    long
    colStride() const { return ind_.cs; }

    size_type
    size() const { return ind_.area(); }

    const MRange&
    ind() const { return ind_; }

    bool
    contiguous() const { return isContiguous(ind_); }

    bool
    transposed() const { return isTransposed(ind_); }

    explicit operator bool() const { return bool(pdata_); }

    pointer
    data() const { return pdata_; }

    void
    applyTrans() { ind_ = MRange(ind_.cn,ind_.cs,ind_.rn,ind_.rs); }

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
    };

std::ostream&
operator<<(std::ostream& s, CMatRef M);

template<typename... CtrArgs>
MatRef
makeMatRef(Real* p,
           CtrArgs&&... args)
    {
    return MatRef(p,std::forward<CtrArgs>(args)...);
    }

template<typename... CtrArgs>
CMatRef
makeMatRef(const Real* p,
           CtrArgs&&... args)
    {
    return CMatRef(p,std::forward<CtrArgs>(args)...);
    }

template<typename MatT>
auto
makeMatRef(MatT& M)
    {
    return makeMatRef(M.data(),M.ind());
    }

template<typename MatT>
auto
makeMatRef(MatT& M,
           const MRange& ind)
    {
    return makeMatRef(M.data(),ind);
    }

template<typename D>
MatRefT<D>
makeMatRef(MatRefT<D> M) 
    { 
    return M; 
    }

template<typename D>
MatRefT<D>
makeMatRef(MatRefT<D> M,
           const MRange& ind)
    {
    return MatRefT<D>(M.data(),ind);
    }


//Copy data referenced by b to memory referenced by a
void
operator&=(const MatRef& a, CMatRef b);

void
operator*=(const MatRef& v, Real fac);

void
operator/=(const MatRef& v, Real fac);

void
operator+=(const MatRef& a, CMatRef b);

void
operator-=(const MatRef& a, CMatRef b);

// compute matrix multiply (dgemm) A*B
// write result to memory referenced by C
void
mult(CMatRef A, 
     CMatRef B, 
     MatRef  C);

// compute matrix multiply (dgemm) A*B
// add result to memory referenced by C
void
multAdd(CMatRef A, 
        CMatRef B, 
        MatRef  C);

void
mult(CMatRef M,
     CVecRef x,
     VecRef y,
     bool fromleft = false);

Real
norm(CMatRef M);

MatRef
randomize(MatRef M);

class Mat
    {
    public:
    using storage_type = std::vector<Real>;
    using iterator = storage_type::iterator;
    using const_iterator = storage_type::const_iterator;
    using value_type = Real;
    using size_type = storage_type::size_type;
    public:
    MRange ind_;
    storage_type data_;
    public:

    Mat() { }

    Mat(long nrows,
        long ncols) 
        : 
        ind_(nrows,ncols),
        data_(nrows*ncols,0) 
        { }

    Mat(const Mat& other) { assignFromMat(other); }

    Mat(Mat&& other) { moveFromMat(std::move(other)); }

    explicit
    Mat(CMatRef ref) { assignFromRef(ref); }

    Mat&
    operator=(const Mat& other) { assignFromMat(other); return *this; }
    Mat& 
    operator=(Mat&& other) { moveFromMat(std::move(other)); return *this; }
    Mat&
    operator=(CMatRef ref) { assignFromRef(ref); return *this; }

    operator CMatRef() const { return CMatRef(data_.data(),ind_); }

    long
    Nrows() const { return ind_.rn; }
    long
    Ncols() const { return ind_.cn; }
    long
    rowStride() const { return ind_.rs; }
    long
    colStride() const { return ind_.cs; }

    const MRange&
    ind() const { return ind_; }

    bool
    contiguous() const { return true; }

    bool
    transposed() const { return false; }

    Real&
    operator()(long i, long j)
        { 
#ifdef DEBUG
        return data_.at(ind_.index(i,j)); 
#else
        return data_[ind_.index(i,j)]; 
#endif
        }

    Real
    operator()(long i, long j) const
        { 
#ifdef DEBUG
        return data_.at(ind_.index(i,j)); 
#else
        return data_[ind_.index(i,j)]; 
#endif
        }

    Mat&
    operator*=(Real fac);
    Mat&
    operator/=(Real fac);
    Mat&
    operator+=(const Mat& other);
    Mat&
    operator-=(const Mat& other);
    Mat&
    operator+=(CMatRef other);
    Mat&
    operator-=(CMatRef other);

    explicit operator bool() const { return !data_.empty(); }

    Real*
    data() { return data_.data(); }

    const Real*
    data() const { return data_.data(); }

    size_type
    size() const { return data_.size(); }

    void
    clear() 
        { 
        ind_ = MRange();
        data_.clear(); 
        }

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
    assignFromRef(CMatRef other)
        {
        ind_ = MRange(other.Nrows(),other.Ncols());
        //Copy data from other contiguously into data_
        data_ = storage_type(other.cbegin(),other.cend());
        }

    void
    assignFromMat(const Mat& other)
        {
        if(&other == this) return;
        ind_ = other.ind_;
        data_ = other.data_;
        }

    void
    moveFromMat(Mat&& other)
        {
        ind_ = other.ind_;
        data_ = std::move(other.data_);
        other.clear();
        }

    };

Mat inline
operator+(Mat A, const Mat& B) { A += B; return A; }
Mat inline
operator+(Mat A, CMatRef B) { A += B; return A; }
Mat inline
operator+(const Mat& A, Mat&& B)
    {
    Mat res(std::move(B));
    res += A;
    return res;
    }
Mat inline
operator-(Mat A, const Mat& B) { A -= B; return A; }
Mat inline
operator-(Mat A, CMatRef B) { A -= B; return A; }
Mat inline
operator-(const Mat& A, Mat&& B)
    {
    Mat res(std::move(B));
    res *= -1;
    res += A;
    return res;
    }
Mat inline
operator*(Mat A, Real fac) { A *= fac; return A; }
Mat inline
operator*(Real fac, Mat A) { return operator*(A,fac); }
Mat inline
operator/(Mat A, Real fac) { A /= fac; return A; }

void inline
mult(CMatRef A,
     CMatRef B,
     Mat& C)
    { 
    mult(A,B,makeMatRef(C)); 
    }

void inline
mult(CMatRef M,
     CVecRef x,
     Vec& y,
     bool fromleft = false)
    { 
    mult(M,x,makeVecRef(y),fromleft); 
    }

Mat inline
operator*(CMatRef A,
          CMatRef B)
    {
    Mat C(A.Nrows(),B.Ncols());
    mult(A,B,C);
    return C;
    }

Vec inline
operator*(CMatRef A,
          CVecRef v)
    {
    Vec res(A.Nrows());
    mult(A,v,res);
    return res;
    }

Vec inline
operator*(CVecRef v,
          CMatRef A)
    {
    Vec res(A.Ncols());
    bool fromleft = true;
    mult(A,v,res,fromleft);
    return res;
    }

inline Mat&
randomize(Mat& v) { randomize(makeMatRef(v)); return v; }

template<typename... CtrArgs>
Mat
randomMat(CtrArgs&&... args)
    {
    Mat M(std::forward<CtrArgs>(args)...);
    randomize(M);
    return M;
    }


};

#endif
