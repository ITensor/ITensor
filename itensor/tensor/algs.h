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
#ifndef __ITENSOR_MATRIX_ALGS__H_
#define __ITENSOR_MATRIX_ALGS__H_

#include "itensor/tensor/slicemat.h"
#include "itensor/util/args.h"

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
    const Args & args = Args::global() );

  
//
// Compute QR decomposition of MxN A matrix such that 
// norm(A-QR) < epsilon
// where Q is MxM orthogonal (unitary) matrix
// and R is MxN upper triangular.
// If complete = false, instead compute "thin" QR: for M >= N
//   Q is MxN matrix with orthonormal columns: Q^T Q = I
//   R is NxN upper triangular
template<class MatA, 
         class MatQ,
         class MatR>
void
QR(MatA && A,
   MatQ && Q,
   MatR && R,
   const Args & args = Args::global());

//
// Hermitian Matrix exponentiate
// by diagHermitian
//
template<class MatM,
         class ScalarT,
         class = stdx::require<hasMatRange<MatM>>>
Mat<common_type<val_type<MatM>,ScalarT>>
expHermitian(MatM && M,
                   ScalarT t);

template<class MatM,
         class = stdx::require<hasMatRange<MatM>>>
Mat<val_type<MatM>>
expHermitian(MatM && M) { return expHermitian(M,1.); }

//
// expMatrix computes exp(tH),
// where H is a general dense matrix,
// and t can be real or complex.
// Either of the two options can be used:
// 1. the irreducible rational Pade approximation
// to the exponential exp(z) = r(z) = (+/-)(I+2*(q(z)/p(z))),
// combined with scaling and squaring;
// 2. the uniform rational Chebyshev approximation to exp(-x) of type(14,14),
// using which 14-digit accuracy is expected if tH is negative definite,
// but may behave poorly otherwise.
//
template<class MatM,
         class ScalarT,
         class = stdx::require<hasMatRange<MatM>>>
Mat<common_type<val_type<MatM>,ScalarT>>
expMatrix(MatM && M,
          ScalarT t,
          int ideg = 6);

template<class MatM,
         class = stdx::require<hasMatRange<MatM>>>
Mat<val_type<MatM>>
expMatrix(MatM && M) { return expMatrix(M,1.); }

} //namespace itensor

#include "itensor/tensor/algs_impl.h"

#endif
