//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MATRIX_ALGS__H_
#define __ITENSOR_MATRIX_ALGS__H_

#include "matrix.h"

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
//   M == U*D*U.t(); //<-- pseudo code
//
// (Note this is the transpose of the convention
//  for U used in the EigenValues routine of
//  earlier version of ITensor/MatrixRef)
//

void
diagSymmetric(const matrixref& M,
              matrixref& U,
              vecref& d);

void
diagSymmetric(const matrixref& M,
              matrix& U,
              vec& d);

//orthogonalize the first num columns of a matrixref M,
//optionally repeating numpass times to reduce roundoff errors
void 
orthog(const matrixref& M, long num = -1, long numpass = 2);

void
SVD(const matrixref& A,
    matrix& U, 
    vec& D, 
    matrix& V,
    Real thresh = 1E-3);

};

#endif
