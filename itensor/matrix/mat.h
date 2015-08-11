//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MAT__H_
#define __ITENSOR_MAT__H_

#include "itensor/matrix/vec.h"
#include "itensor/matrix/matiter.h"

namespace itensor {

template<typename T>
class Matrix;

template<typename T>
class MatrixRef;

using Mat = Matrix<Real>;
using MatRef = MatrixRef<Real>;
using MatRefc = MatrixRef<const Real>;

//using CMat = Matrix<Complex>;
//using CMatRef = MatrixRef<Complex>;
//using CMatRefc = MatrixRef<const Complex>;

template<typename T>
class MatrixRef
    {
    public:
    using iterator = MatIter<T*>;
    using const_iterator = MatIter<const T*>;
    using value_type = std::remove_const_t<T>;
    using pointer = T*;
    using reference = T&;
    using size_type = long;
    using mat_type = std::conditional_t<std::is_const<T>::value,
                                        const Matrix<value_type>,
                                        Matrix<value_type>>;
    private:
    pointer pdata_ = nullptr;
    MRange range_;
    public:

    MatrixRef() { }

    MatrixRef(pointer pdata, 
              long nrows,
              long ncols)
        :
        pdata_(pdata),
        range_(nrows,ncols)
        { }


    MatrixRef(pointer pdata, 
              MRange const& ind)
        :
        pdata_(pdata),
        range_(ind)
        { }

    //MatrixRef(pointer pdata, 
    //          long offset,
    //          long nrows,
    //          long ncols,
    //          bool trans)
    //    : MatrixRef(pdata+offset,nrows,ncols,trans)
    //    { }

    MatrixRef(pointer pdata, 
              long offset,
              MRange const& ind)
        : MatrixRef(pdata+offset,ind)
        { }

    MatrixRef(mat_type& M) { pointTo(M); }

    MatrixRef&
    operator=(mat_type& M) { pointTo(M); return *this; }

    MatrixRef(mat_type&& M) = delete;

    MatrixRef&
    operator=(mat_type&& M) = delete;

    operator MatrixRef<const T>() const { return MatrixRef<const T>(pdata_,range_); }

    explicit operator bool() const { return bool(pdata_); }

    long
    Nrows() const { return range_.rn; }
    long
    Ncols() const { return range_.cn; }
    long
    rowStride() const { return range_.rs; }
    long
    colStride() const { return range_.cs; }

    size_type
    size() const { return range_.area(); }

    MRange const&
    range() const { return range_; }

    bool
    contiguous() const { return isContiguous(range_); }

    bool
    transposed() const { return isTransposed(range_); }

    pointer
    data() const { return pdata_; }

    reference
    operator()(long i, long j) const { return pdata_[range_.index(i,j)]; }

    iterator 
    begin() const { return iterator(pdata_,range_); }

    iterator 
    end() const { return iterator(range_); }

    const_iterator 
    cbegin() const { return const_iterator(pdata_,range_); }

    const_iterator 
    cend() const { return const_iterator(range_); }

    void
    clear()
        {
        pdata_ = nullptr;
        range_ = MRange();
        }

    void
    applyTrans() { range_ = MRange(range_.cn,range_.cs,range_.rn,range_.rs); }

    private:
    void
    pointTo(mat_type& M);
    };


void
operator*=(const MatRef& v, Real fac);

void
operator/=(const MatRef& v, Real fac);

void
operator+=(const MatRef& a, MatRefc b);

void
operator-=(const MatRef& a, MatRefc b);

//Copy data referenced by b to memory referenced by a
void
operator&=(const MatRef& a, MatRefc b);

void
call_dgemm(MatRefc A, 
           MatRefc B, 
           MatRef  C,
           Real alpha,
           Real beta);

void
mult(MatRefc A,
     MatRefc B,
     MatRef  C);


// compute matrix multiply (dgemm) A*B
// add result to memory referenced by C
void
multAdd(MatRefc A, 
        MatRefc B, 
        MatRef  C);

void
mult(MatRefc M,
     VecRefc x,
     VecRef y,
     bool fromleft = false);


template<typename T>
class Matrix
    {
    public:
    using storage_type = std::vector<T>;
    using iterator = typename storage_type::iterator;
    using const_iterator = typename storage_type::const_iterator;
    using value_type = std::remove_const_t<T>;
    using pointer = std::add_pointer<T>;
    using reference = std::add_lvalue_reference_t<T>;
    using size_type = long;
    public:
    MRange range_;
    storage_type data_;
    public:

    Matrix() { }

    Matrix(long nrows,
           long ncols) 
        : 
        range_(nrows,ncols),
        data_(nrows*ncols,0) 
        { }

    Matrix(long nrows,
           long ncols,
           value_type val)
        : 
        range_(nrows,ncols),
        data_(nrows*ncols,val) 
        { }

    Matrix(const Matrix& other) { assignFromMat(other); }

    Matrix(Matrix&& other) { moveFromMat(std::move(other)); }

    explicit
    Matrix(const MatrixRef<value_type>& ref) { assignFromRef(ref); }

    explicit
    Matrix(const MatrixRef<const value_type>& ref) { assignFromRef(ref); }

    Matrix&
    operator=(const Matrix& other) { assignFromMat(other); return *this; }

    Matrix& 
    operator=(Matrix&& other) { moveFromMat(std::move(other)); return *this; }

    Matrix&
    operator=(const MatrixRef<value_type>& ref) { assignFromRef(ref); return *this; }

    Matrix&
    operator=(const MatrixRef<const value_type>& ref) { assignFromRef(ref); return *this; }

    explicit operator bool() const { return !data_.empty(); }

    size_type
    size() const { return data_.size(); }

    long
    Nrows() const { return range_.rn; }
    long
    Ncols() const { return range_.cn; }
    long
    rowStride() const { return range_.rs; }
    long
    colStride() const { return range_.cs; }

    const MRange&
    range() const { return range_; }

    bool
    contiguous() const { return true; }

    bool
    transposed() const { return false; }

    reference
    operator()(long i, long j)
        { 
#ifdef DEBUG
        return data_.at(range_.index(i,j)); 
#else
        return data_[range_.index(i,j)]; 
#endif
        }

    value_type
    operator()(long i, long j) const
        { 
#ifdef DEBUG
        return data_.at(range_.index(i,j)); 
#else
        return data_[range_.index(i,j)]; 
#endif
        }

    Matrix&
    operator*=(Real fac);

    Matrix&
    operator/=(Real fac);

    Matrix&
    operator+=(const Matrix& other);

    Matrix&
    operator-=(const Matrix& other);

    Matrix&
    operator+=(MatrixRef<const T> other);

    Matrix&
    operator-=(MatrixRef<const T> other);

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

    T*
    data() { return data_.data(); }

    const T*
    data() const { return data_.data(); }

    storage_type&
    store() { return data_; }

    const storage_type&
    store() const { return data_; }

    //Essentially no cost when resizing downward
    void
    resize(long nrows,
           long ncols)
        {
        range_ = MRange(nrows,ncols);
        data_.resize(range_.area(),0);
        }

    //Reducing number of columns does not affect
    //remaining data (column major storage)
    void
    reduceColsTo(long newcols)
        {
#ifdef DEBUG
        if(newcols > Ncols()) throw std::runtime_error("newcols must be less than current Ncols()");
#endif
        range_ = MRange(Nrows(),newcols);
        data_.resize(range_.area());
        }

    void
    clear() 
        { 
        range_ = MRange();
        data_.clear(); 
        }

    private:

    void
    assignFromRef(const MatrixRef<const value_type>& other)
        {
        range_ = MRange(other.Nrows(),other.Ncols());
        //Copy data from other contiguously into data_
        data_ = storage_type(other.cbegin(),other.cend());
        }

    void
    assignFromMat(const Matrix& other)
        {
        if(&other == this) return;
        range_ = other.range_;
        data_ = other.data_;
        }

    void
    moveFromMat(Matrix&& other)
        {
        range_ = other.range_;
        data_ = std::move(other.data_);
        other.clear();
        }
    };

template<typename T>
void MatrixRef<T>::
pointTo(mat_type& M)
    {
    pdata_ = M.data();
    range_ = M.range();
    }

//
// makeRef functions
//

template<typename T>
auto
makeRef(MatrixRef<T>& M) { return M; }

template<typename T>
auto
makeRef(MatrixRef<const T>& M) { return M; }

template<typename T>
auto
makeRef(Matrix<T>& M) { return MatrixRef<T>(M); }

template<typename T>
auto
makeRef(const Matrix<T>& M) { return MatrixRef<const T>(M); }

template<typename T, typename Arg, typename... Rest>
auto
makeRef(const MatrixRef<T>& M, Arg&& arg, Rest&&... args) 
    { return MatrixRef<T>(M.data(),std::forward<Arg>(arg),std::forward<Rest>(args)...); }

template<typename T, typename Arg, typename... Rest>
auto
makeRef(Matrix<T>& M, Arg&& arg, Rest&&... args) 
    { return MatrixRef<T>(M.data(),std::forward<Arg>(arg),std::forward<Rest>(args)...); }

template<typename T, typename Arg, typename... Rest>
auto
makeRef(const Matrix<T>& M, Arg&& arg, Rest&&... args) 
    { return MatrixRef<const T>(M.data(),std::forward<Arg>(arg),std::forward<Rest>(args)...); }


//This version of makeRef intended to fail,
//forbids explicitly making MatrixRef's to temporaries
template<typename T, typename... Rest>
auto
makeRef(Matrix<T>&& M, Rest&&... args) { return MatrixRef<T>(std::move(M)); }


//
// makeVecRef
//

template<typename T>
auto
makeVecRef(Matrix<T>& M) { return VectorRef<T>(M.data(),M.size()); }

template<typename T>
auto
makeVecRef(const Matrix<T>& M) { return VectorRef<const T>(M.data(),M.size()); }

template<typename T, typename Arg, typename... Rest>
auto
makeVecRef(Matrix<T>& M, Arg&& arg, Rest&&... args) 
    { return VectorRef<T>(M.data(),std::forward<Arg>(arg),std::forward<Rest>(args)...); }

template<typename T, typename Arg, typename... Rest>
auto
makeVecRef(const Matrix<T>& M, Arg&& arg, Rest&&... args) 
    { return VectorRef<const T>(M.data(),std::forward<Arg>(arg),std::forward<Rest>(args)...); }

template<typename T, typename Arg, typename... Rest>
auto
makeVecRef(const MatrixRef<T>& M, Arg&& arg, Rest&&... args) 
    { return VectorRef<T>(M.data(),std::forward<Arg>(arg),std::forward<Rest>(args)...); }

//This version of makeVecRef intended to fail,
//forbids explicitly making VectorRef's to temporaries
template<typename T, typename... Rest>
auto
makeVecRef(Matrix<T>&& M, Rest&&... args) { return VectorRef<T>(std::move(M)); }


//
// Finish defining Matrix operators
//

template<typename T>
Matrix<T>& Matrix<T>::
operator*=(Real fac) { makeRef(*this) *= fac; return *this; }

template<typename T>
Matrix<T>& Matrix<T>::
operator/=(Real fac) { makeRef(*this) /= fac; return *this; }

template<typename T>
Matrix<T>& Matrix<T>::
operator+=(const Matrix<T>& other) { makeRef(*this) += makeRef(other); return *this; }

template<typename T>
Matrix<T>& Matrix<T>::
operator-=(const Matrix<T>& other) { makeRef(*this) -= makeRef(other); return *this; }

template<typename T>
Matrix<T>& Matrix<T>::
operator+=(MatrixRef<const T> other) {  makeRef(*this) += other; return *this; }

template<typename T>
Matrix<T>& Matrix<T>::
operator-=(MatrixRef<const T> other) { makeRef(*this) -= other; return *this; }

template<typename T>
Matrix<T>
operator*(Matrix<T> A, Real fac) { A *= fac; return A; }

template<typename T>
Matrix<T>
operator*(Real fac, Matrix<T> A) { A *= fac; return A; }

template<typename T>
Matrix<T>
operator/(Matrix<T> A, Real fac) { A /= fac; return A; }

Mat inline
operator+(MatRefc A, MatRefc B)
    { 
    Mat res(A);
    res += B; 
    return res; 
    }

Mat inline
operator+(MatRefc A, Mat&& B) 
    { 
    Mat res(std::move(B)); 
    res += A; 
    return res; 
    }

Mat inline
operator+(Mat&& A, MatRefc B) 
    { 
    Mat res(std::move(A)); 
    res += B; 
    return res; 
    }

Mat inline
operator-(MatRefc A, MatRefc B)
    { 
    Mat res(A);
    res -= B; 
    return res; 
    }

Mat inline
operator-(MatRefc A, Mat&& B) 
    { 
    Mat res(std::move(B)); 
    res *= -1;
    res += A; 
    return res; 
    }

Mat inline
operator-(Mat&& A, MatRefc B) 
    { 
    Mat res(std::move(A)); 
    res -= B; 
    return res; 
    }

Mat inline
matrixMult(MatRefc A,
           MatRefc B)
    {
    Mat C(A.Nrows(),B.Ncols());
    call_dgemm(A,B,C,1.,0.);
    return C;
    }

Mat inline
operator*(MatRefc A, const Mat& B) { return matrixMult(A,B); }

Mat inline
operator*(const Mat& A, MatRefc B) { return matrixMult(A,B); }

Mat inline
operator*(const Mat& A, const Mat& B) { return matrixMult(A,B); }

Mat inline
operator*(MatRefc A, MatRefc B) { return matrixMult(A,B); }

inline Mat&
operator*=(Mat& A, MatRefc B) { A = matrixMult(A,B); return A; }

Vec inline
operator*(MatRefc A,
          VecRefc v)
    {
    Vec res(A.Nrows());
    mult(A,v,res);
    return res;
    }

Vec inline
operator*(VecRefc v,
          MatRefc A)
    {
    Vec res(A.Ncols());
    bool fromleft = true;
    mult(A,v,res,fromleft);
    return res;
    }

Real
norm(MatRefc M);

Real inline
norm(const Mat& M) { return norm(makeRef(M)); }

//
// These versions of op-assign to MatRef can
// work safely for temporary Mat's since
// const references extend lifetime of rvalues
//

void inline
operator&=(MatRef a, Mat const& b) { a &= makeRef(b); }

void inline
operator+=(MatRef a, Mat const& b) { a += makeRef(b); }

void inline
operator-=(MatRef a, Mat const& b) { a -= makeRef(b); }

void
randomize(MatRef M);

std::ostream&
operator<<(std::ostream& s, MatRefc M);

inline std::ostream&
operator<<(std::ostream& s, Mat const& M) { return operator<<(s,makeRef(M)); }

template<typename... CtrArgs>
Mat
randomMat(CtrArgs&&... args)
    {
    Mat M(std::forward<CtrArgs>(args)...);
    randomize(M);
    return M;
    }

} //namespace itensor

#endif
