//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MATRIX_ALGS__H_
#define __ITENSOR_MATRIX_ALGS__H_

#include "itensor/tensor/slicemat.h"

namespace itensor {

//
// diagSymmetric diagonalizes a (real)
// symmetric matrix M and return U, d
// such that:
// o Elements of d are the eigenvalues in 
//   decreasing order.
// o Columns of U are the corresponding eigenvectors.
//
// In other words, the following holds:
//
//   diagSymmetric(M,U,d);
//   auto D = Matrix(nrows(M),ncols(M));
//   diagonal(D) &= d;
//   M == U*D*transpose(U); //<-- pseudo code
//
//

void
diagSymmetric(MatrixRefc const& M,
              Matrix          & U,
              Vector          & d);

void
diagSymmetric(MatrixRefc const& M,
              Matrix          & U,
              VectorRef  const& d);

void
diagSymmetric(MatrixRefc const& M,
              MatrixRef  const& U,
              VectorRef  const& d);

//
// diagHermitian diagonalizes a 
// hermitian matrix with real part Mre (symmetric)
// and imaginary part Mim (antisymmetric)
// such that:
// o Elements of d are the (real) eigenvalues in 
//   decreasing order.
// o Columns of Ure and Uim are the real and imag
//   parts of the corresponding eigenvectors
//
// In other words, the following holds:
//
//   diagHermitian(Mre,Mim,Ure,Uim,d);
//   auto D = Matrix(nrows(Mre),ncols(Mre));
//   diagonal(D) &= d;
//   Mre == Ure*D*transpose(Ure) + Uim*D*transpose(Uim); //<-- pseudo code
//   Mim == Uim*D*transpose(Ure) - Ure*D*transpose(Uim); //<-- pseudo code
//

void
diagHermitian(MatrixRefc const& Mre,
              MatrixRefc const& Mim,
              MatrixRef  const& Ure,
              MatrixRef  const& Uim,
              VectorRef  const& d);

void
diagHermitian(MatrixRefc const& Mre,
              MatrixRefc const& Mim,
              Matrix          & Ure,
              Matrix          & Uim,
              VectorRef  const& d);

void
diagHermitian(MatrixRefc const& Mre,
              MatrixRefc const& Mim,
              Matrix          & Ure,
              Matrix          & Uim,
              Vector          & d);

//orthogonalize the first num columns of a matrixref M,
//optionally repeating numpass times to reduce roundoff errors
//if num == 0, orthogonalizes all columns
void 
orthog(MatrixRef M, size_t num = 0, size_t numpass = 2);

//Compute U,D,V such that 
//norm(A-U*DD*transpose(V)) < epsilon
//where DD is matrix with D on diagonal
void
SVDRef(MatrixRefc const& M,
       MatrixRef  const& U, 
       VectorRef  const& D, 
       MatrixRef  const& V,
       Real thresh = 1E-3);

//Compute U,D,V such that 
//norm(A-U*D*transpose(V)) < epsilon
//where DD is matrix with D on diagonal
void
SVD(MatrixRefc const& A,
    Matrix & U, 
    Vector & D, 
    Matrix & V,
    Real thresh = 1E-3);

void 
checksvd(MatrixRefc const& A, 
         MatrixRefc const& U, 
         VectorRefc const& D, 
         MatrixRefc const& V);

} //namespace itensor

#endif
