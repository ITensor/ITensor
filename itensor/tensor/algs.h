//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
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

template<class MatA, class MatB, class MatX,
         class = stdx::require<
         hasMatRange<MatA>,
         hasMatRange<MatB>,
         hasMatRange<MatX>
         >>
void
linSystem(MatA && A,
          MatB && B,
          MatX && X,
          Args const& args);

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


namespace detail {
  //Helper for linSystem
  template<typename Iter, typename V>
  void
  copyElts(Iter m,
          MatRef<V> const& U)
  {
    auto ue = U.data()+U.size();
    for(auto u = U.data(); u != ue; ++u,++m)
      {
      *u = (*m);
      }
  }

  int
  hermitianDiag(int N, Real *Udata, Real *ddata);
  int
  hermitianDiag(int N, Cplx *Udata,Real *ddata);

  int
  dtl_linSystem(int N, int Nb, Real const * Adata, Real * Bdata, Args const& args);
  int
  dtl_linSystem(int N, int Nb, Cplx const * Adata, Cplx * Bdata, Args const& args);
} //namespace detail

template<class MatA, 
         class MatB,
         class MatX,
         class>
void
linSystem(MatA && A,
          MatB && B,
          MatX && X,
          Args const& args)
{
  using Aval = typename stdx::decay_t<MatA>::value_type;
  using Bval = typename stdx::decay_t<MatB>::value_type;
  static_assert((isReal<Aval>() && isReal<Bval>()) || (isCplx<Aval>() && isCplx<Bval>()),
              "A and B must be both real or both complex in linSystem");
  
  auto N = ncols(A);
  if(N < 1) throw std::runtime_error("linSystem: 0 dimensional matrix");
  if(N != nrows(A))
    {
    printfln("A is %dx%d",nrows(A),ncols(A));
    throw std::runtime_error("linSystem: Input Matrix must be square");
    }
  
  // copy B into X
  auto Nb = ncols(B);
  resize(X,nrows(B),Nb);
  
  detail::copyElts(B.data(),makeRef(X));

  auto info = detail::dtl_linSystem(N,Nb,A.data(),X.data(),args);
  if(info != 0) 
    {
      std::cout<<"[linSystem] Error info= "<< info << std::endl;
      auto arg_method = args.getString("method","LU");
      if( (info > 0) && (arg_method=="CHOLESKY") ) {
        printfln("Minor %d not positive",info);

        std::vector<double> d(N);
        if(isCplx<Aval>()) {
          auto pA = reinterpret_cast<Cplx const*>(A.data());

          std::vector<Cplx> cpA;
          cpA.resize(N*N);
          std::copy(pA,pA+N*N,cpA.data());

          auto diagInfo = detail::hermitianDiag(N,cpA.data(), d.data());
        } else {
          auto pA = reinterpret_cast<double const*>(A.data());

          std::vector<double> cpA;
          cpA.resize(N*N);
          std::copy(pA,pA+N*N,cpA.data());

          auto diagInfo = detail::hermitianDiag(N,cpA.data(), d.data());          
        }
        std::cout<<"Spectrum: ";
        for (auto const& elem : d) { std::cout<< elem <<" "; }
        std::cout<<std::endl;
        throw std::runtime_error("Error condition in linSystem");
      } else {
        println("A = \n",A);
        throw std::runtime_error("Error condition in linSystem");
      }
    }
}


} //namespace itensor

#include "itensor/tensor/algs.ih"

#endif
