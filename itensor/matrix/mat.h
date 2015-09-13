//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MAT__H_
#define __ITENSOR_MAT__H_

#include "itensor/matrix/matrange.h"

namespace itensor {

using MatrixRef = TenRef<MatRange>;
using MatrixRefc = TenRefc<MatRange>;
using Matrix = Ten<MatRange>;

template<typename Mat_>
auto
nrows(Mat_ const& M) { return M.range().rn; }

template<typename Mat_>
auto
ncols(Mat_ const& M) { return M.range().cn; }

template<typename Mat_>
auto
rowStride(Mat_ const& M) { return M.range().rs; }

template<typename Mat_>
auto
colStride(Mat_ const& M) { return M.range().cs; }

template<typename Mat_>
auto
isTransposed(Mat_ const& M) { return isTransposed(M.range()); }


void
operator*=(MatrixRef const& M, Real fac);

void
operator/=(MatrixRef const& M, Real fac);

void
operator+=(MatrixRef const& A, MatrixRefc const& B);

void
operator-=(MatrixRef const& A, MatrixRefc const& B);

//Copy data referenced by B to memory referenced by A
void
operator&=(MatrixRef const& A, MatrixRefc const& B);

//Copy data of B to memory referenced by A
void inline
operator&=(MatrixRef const& A, Matrix const& B) { A &= makeRefc(B); }

void
call_dgemm(MatrixRefc A, 
           MatrixRefc B, 
           MatrixRef  C,
           Real alpha,
           Real beta);

void
mult(MatrixRefc A,
     MatrixRefc B,
     MatrixRef  C);


// compute matrix multiply (dgemm) A*B
// add result to memory referenced by C
void
multAdd(MatrixRefc A, 
        MatrixRefc B, 
        MatrixRef  C);

void
mult(MatrixRefc M,
     VectorRefc x,
     VectorRef y,
     bool fromleft = false);

//y = y+M*x
void
multAdd(MatrixRefc M,
        VectorRefc x,
        VectorRef y,
        bool fromleft = false);

//y = y-M*x
void
multSub(MatrixRefc M,
        VectorRefc x,
        VectorRef y,
        bool fromleft = false);

//Reducing number of columns does not affect
//remaining data (column major storage)
void
reduceCols(Matrix & M, size_t new_ncols);

void
resize(Matrix & M, size_t nrows, size_t ncols);


Matrix 
operator*(MatrixRefc const& A, Real fac);

Matrix 
operator*(Real fac, MatrixRefc const& A);

Matrix 
operator*(Matrix && A, Real fac);

Matrix 
operator*(Real fac, Matrix && A);

Matrix 
operator/(MatrixRefc const& A, Real fac);

Matrix 
operator/(Matrix && A, Real fac);

Matrix 
operator+(MatrixRefc const& A, MatrixRefc const& B);

Matrix 
operator+(MatrixRefc const& A, Matrix && B);

Matrix 
operator+(Matrix && A, MatrixRefc const& B);

Matrix 
operator-(MatrixRefc const& A, MatrixRefc const& B);

Matrix 
operator-(MatrixRefc const& A, Matrix && B);

Matrix 
operator-(Matrix && A, MatrixRefc const& B);

Matrix 
matrixMult(MatrixRefc const& A,
           MatrixRefc const& B);

Matrix inline
operator*(MatrixRefc const& A, MatrixRefc const& B) { return matrixMult(A,B); }

Matrix inline
operator*(Matrix const& A, MatrixRefc const& B) { return matrixMult(A,B); }

Matrix inline
operator*(MatrixRefc const& A, Matrix const& B) { return matrixMult(A,B); }

Matrix inline
operator*(Matrix const& A, Matrix const& B) { return matrixMult(A,B); }

inline Matrix&
operator*=(Matrix & A, MatrixRefc const& B) { A = matrixMult(A,B); return A; }

Vector
operator*(MatrixRefc const& A,
          VectorRefc const& v);

Vector
operator*(VectorRefc const& v,
          MatrixRefc const& A);


void
randomize(MatrixRef const& M);
void
randomize(Matrix & M);

template<>
std::ostream&
operator<<(std::ostream& s, MatrixRefc const& M);

template<> inline
std::ostream&
operator<<(std::ostream& s, MatrixRef const& M) { return s << makeRefc(M); }

inline std::ostream&
operator<<(std::ostream& s, Matrix const& M) { return s << makeRefc(M); }

template<typename... CtrArgs>
Matrix
randomMat(CtrArgs&&... args)
    {
    Matrix M(std::forward<CtrArgs>(args)...);
    randomize(M);
    return M;
    }

} //namespace itensor

#endif
