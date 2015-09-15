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
//   matrix D = matrix(M.Nrows(),M.Ncols());
//   diagonal(D) = d;
//   M == U*D*transpose(U); //<-- pseudo code
//
// (Note this is the transpose of the convention
//  for U used in the EigenValues routine of
//  earlier version of ITensor/MatrixRef)
//

//void
//diagSymmetric(const MatrixRefc& M,
//              MatrixRef&  U,
//              VectorRef&  d);

void
diagSymmetric(MatrixRefc M,
              Matrix& U,
              Vector& d);

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
