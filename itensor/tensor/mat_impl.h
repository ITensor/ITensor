//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MAT_IMPL_H__
#define __ITENSOR_MAT_IMPL_H__

namespace itensor {

template<typename V>
void
reduceCols(Mat<V> & M, size_t new_ncols)
    {
#ifdef DEBUG
    if(new_ncols > ncols(M)) throw std::runtime_error("new ncols > old ncols in reduceCols");
#endif
    M.resize(MatRange(nrows(M),new_ncols));
    }

template<typename V>
void
resize(Mat<V> & M, size_t nrows, size_t ncols)
    {
    M.resize(MatRange(nrows,ncols));
    }

template<typename MatA,
         class = stdx::require<hasMatRange<MatA>> >
Mat<val_type<MatA>> 
operator*(MatA const& A, Real fac)
    { 
    Mat<val_type<MatA>> res(A);
    res *= fac; 
    return res; 
    }

template<typename MatA,
         class = stdx::require<hasMatRange<MatA>> >
Mat<val_type<MatA>> 
operator*(Real fac, MatA const& A)
    { 
    Mat<val_type<MatA>> res(A);
    res *= fac; 
    return res; 
    }

template<typename V>
Mat<V> 
operator*(Mat<V> && A, Real fac)
    { 
    Mat<V> res(std::move(A));
    res *= fac; 
    return res; 
    }

template<typename V>
Mat<V> 
operator*(Real fac, Mat<V> && A)
    { 
    Mat<V> res(std::move(A));
    res *= fac; 
    return res; 
    }

template<typename MatA,
         class = stdx::require<hasMatRange<MatA>> >
Mat<val_type<MatA>> 
operator/(MatA const& A, Real fac)
    { 
    Mat<val_type<MatA>> res(A);
    res /= fac; 
    return res; 
    }

template<typename V>
Mat<V> 
operator/(Mat<V> && A, Real fac)
    { 
    Mat<V> res(std::move(A));
    res /= fac; 
    return res; 
    }

template<typename MatA,
         typename MatB,
         class>
Mat<common_type<MatA,MatB>> 
mult(MatA const& A,
     MatB const& B)
    {
    Mat<common_type<MatA,MatB>> C(nrows(A),ncols(B));
    gemm(A,B,makeRef(C),1.,0.);
    return C;
    }

Vector inline
operator*(MatrixRefc const& A,
          VectorRefc const& b)
    {
    Vector res(nrows(A));
    mult(A,b,makeRef(res));
    return res;
    }

Vector inline
operator*(VectorRefc const& a,
          MatrixRefc const& B)
    {
    Vector res(ncols(B));
    bool fromleft = true;
    mult(B,a,makeRef(res),fromleft);
    return res;
    }

template<typename... CtrArgs>
Matrix
randomMat(CtrArgs&&... args)
    {
    auto M = Matrix(std::forward<CtrArgs>(args)...);
    randomize(M);
    return M;
    }

template<typename... CtrArgs>
CMatrix
randomMatC(CtrArgs&&... args)
    {
    auto M = CMatrix(std::forward<CtrArgs>(args)...);
    randomize(M);
    return M;
    }

} //namespace itensor

#endif
