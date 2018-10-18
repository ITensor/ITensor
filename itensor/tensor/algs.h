//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MATRIX_ALGS__H_
#define __ITENSOR_MATRIX_ALGS__H_

#include "itensor/tensor/slicemat.h"

namespace itensor {

static const Real SVD_THRESH = 1E-5;

//
// diagHermitian diagonalizes a
// Hermitian (and/or real symmetric) matrix M 
// and return U, d
// such that:
// o Elements of d are the eigenvalues in 
//   decreasing order.
// o Columns of U are the corresponding eigenvectors.
//
// In other words, the following holds:
//
//   diagHermitian(M,U,d);
//   auto D = Matrix(nrows(M),ncols(M));
//   diagonal(D) &= d;
//   //Cplx case:
//   norm(M-U*D*conj(transpose(U))); //is small
//   //Real case:
//   norm(M-U*D*transpose(U)); //is small
//
template<class MatM, class MatU,class Vecd,
         class = stdx::require<
         hasMatRange<MatM>,
         hasMatRange<MatU>,
         hasVecRange<Vecd>
         >>
void
diagHermitian(MatM && M,
              MatU && U,
              Vecd && d);

// compute eigenvalues
// and right eigenvectors
template<class MatM, class MatV,class Vecd,
         class = stdx::require<
         hasMatRange<MatM>,
         hasMatRange<MatV>,
         hasVecRange<Vecd>
         >>
void
eigen(MatM && M,
      MatV && Vr,
      MatV && Vi,
      Vecd && dr,
      Vecd && di);

// full eigenvalue decomp
// M = L*D*R
template<class MatM, class MatV,class Vecd,
         class = stdx::require<
         hasMatRange<MatM>,
         hasMatRange<MatV>,
         hasVecRange<Vecd>
         >>
void
eigDecomp(MatM && M,
          MatV && Lr,
          MatV && Li,
          Vecd && dr,
          Vecd && di,
          MatV && Rr,
          MatV && Ri);

//orthogonalize the columns of a matrixref M, optionally
//repeating numpass times to reduce roundoff errors
template<typename V>
void 
orthog(MatRef<V> M, size_t numpass = 2);

template<typename Mat_,class=stdx::require<hasMatRange<Mat_>>>
void 
orthog(Mat_ && M, size_t numpass = 2) { orthog(makeRef(M),numpass); }


//
// Compute U,D,V such that 
// norm(A-U*DD*conj(transpose(V))) < epsilon
// where DD is matrix with D on diagonal
// diagonal(DD) &= D;
// (for Real case can leave out conj of V)
//
template<class MatM, class MatU,class VecD,class MatV,
         class = stdx::require<
         hasMatRange<MatM>,
         hasMatRange<MatU>,
         hasVecRange<VecD>,
         hasMatRange<MatV>
         >>
void
SVD(MatM && M,
    MatU && U, 
    VecD && D, 
    MatV && V,
    Real thresh = SVD_THRESH);


} //namespace itensor

#include "itensor/tensor/algs_impl.h"

#endif
