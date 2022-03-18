//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
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
Mat<common_type<val_type<MatA>,Cplx>>
operator*(MatA const& A, Cplx fac)
    {
    Mat<common_type<val_type<MatA>,Cplx>> res(A);
    res *= fac;
    return res;
    }

template<typename MatA,
         class = stdx::require<hasMatRange<MatA>> >
Mat<common_type<val_type<MatA>,Cplx>>
operator*(Cplx fac, MatA const& A) { return A*fac; }

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
operator*(Mat<V> && A, Cplx fac)
    {
    Mat<V> res(std::move(A));
    res *= fac;
    return res;
    }

inline Mat<Cplx>
operator*(Mat<Real> const& A, Cplx fac)
    {
    auto nr = nrows(A);
    auto nc = ncols(A);
    Mat<Cplx> res(nr,nc);
    auto resend = res.data()+res.size();
    auto Adata = A.data();
    for(auto r = res.data(); r != resend; ++r, ++Adata)
        *r = fac * (*Adata);
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

template<typename V>
Mat<common_type<V,Cplx>>
operator*(Cplx fac, Mat<V> const& A) { return A*fac; }

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

CVector inline
operator*(CMatrixRefc const& A,
          CVectorRefc const& b)
    {
    CVector res(nrows(A));
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

CVector inline
operator*(CVectorRefc const& a,
          CMatrixRefc const& B)
    {
    CVector res(ncols(B));
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

template<typename... CtrArgs>
Matrix
randn(CtrArgs&&... args)
    {
    auto M = Matrix(std::forward<CtrArgs>(args)...);
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> normal{0,1};
    for(auto& el : M) el = normal(gen);
    return M;
    }

template<typename... CtrArgs>
CMatrix
randnC(CtrArgs&&... args)
    {
    auto M = CMatrix(std::forward<CtrArgs>(args)...);
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> normal{0,1};
    for(auto& el : M) el = std::complex(normal(gen),normal(gen));
    return M;
    }

} //namespace itensor

#endif
