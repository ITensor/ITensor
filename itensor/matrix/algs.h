//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MATRIX_ALGS__H_
#define __ITENSOR_MATRIX_ALGS__H_

#include "slicemat.h"

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

//void
//diagSymmetric(const MatRefc& M,
//              MatRef&  U,
//              VecRef&  d);

void
diagSymmetric(MatRefc M,
              Mat& U,
              Vec& d);

//orthogonalize the first num columns of a matrixref M,
//optionally repeating numpass times to reduce roundoff errors
void 
orthog(MatRef M, long num = -1, long numpass = 2);

//Compute U,D,V such that 
//A == U*D*transpose(V), schematically
void
SVDRefs(const MatRefc& A,
        const MatRef&  U, 
        const VecRef&  D, 
        const MatRef&  V,
        Real thresh = 1E-3,
        Real northpass = 2);

//Compute U,D,V such that 
//A == U*D*transpose(V), schematically
void
SVD(const MatRefc& A,
    Mat& U, 
    Vec& D, 
    Mat& V,
    Real thresh = 1E-3,
    Real northpass = 2);

}

#endif
