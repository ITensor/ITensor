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
#ifndef __ITENSOR_MATRIX_ALGS_IMPL_H_
#define __ITENSOR_MATRIX_ALGS_IMPL_H_

namespace itensor {

namespace detail {

//Helper for diagHermitian
template<typename Iter, typename V>
void
copyNegElts(Iter m,
            MatRef<V> const& U)
    {
    auto ue = U.data()+U.size();
    for(auto u = U.data(); u != ue; ++u,++m)
        {
        *u = -(*m);
        }
    }

//Helper for diagHermitian
template<typename Iter>
void
copyNegElts(Iter mre,
            Iter mim,
            std::vector<Cplx> & Mc)
    {
    for(auto& z : Mc)
        {
        realRef(z) = -(*mre);
        imagRef(z) = -(*mim);
        ++mre;
        ++mim;
        }
    }
  
  //Helpers for real/complex SVDs
  inline Real
  conjIfCplx(Real n) 
  {
    return n;
      }

   inline Cplx
  conjIfCplx(Cplx n)
  {
    return std::conj(n);
      }
    

  int
  hermitianDiag(int N, Real *Udata, Real *ddata);
  int
  hermitianDiag(int N, Cplx *Udata,Real *ddata);

  int
  QR(int M, int N, int Rrows, Real *Qdata, Real *Rdata);
  int
  QR(int M, int N, int Rrows, Cplx *Qdata,Cplx *Rdata);

  int
  SVD_gesdd(int M, int N, Cplx * Adata, Cplx * Udata, Real * Ddata, Cplx * Vdata);

  int
  SVD_gesdd(int M, int N, Real * Adata, Real * Udata, Real * Ddata, Real * Vdata);

  int
  SVD_gesvd(int M, int N, Cplx * Adata, Cplx * Udata, Real * Ddata, Cplx * Vdata);

  int
  SVD_gesvd(int M, int N, Real * Adata, Real * Udata, Real * Ddata, Real * Vdata);
    

} //namespace detail

template<class MatM, 
         class MatU,
         class Vecd,
         class>
void
diagHermitian(MatM && M,
              MatU && U,
              Vecd && d)
    {
    using Mval = typename stdx::decay_t<MatM>::value_type;
    using Uval = typename stdx::decay_t<MatU>::value_type;
    static_assert((isReal<Mval>() && isReal<Uval>()) || (isCplx<Mval>() && isCplx<Uval>()),
                  "M and U must be both real or both complex in diagHermitian");
    auto N = ncols(M);
    if(N < 1) throw std::runtime_error("diagHermitian: 0 dimensional matrix");
    if(N != nrows(M))
        {
        printfln("M is %dx%d",nrows(M),ncols(M));
        throw std::runtime_error("diagHermitian: Input Matrix must be square");
        }

    resize(U,nrows(M),ncols(M));
    resize(d,nrows(M));

#ifdef DEBUG
    if(!isContiguous(U))
        throw std::runtime_error("diagHermitian: U must be contiguous");
    if(!isContiguous(d))
        throw std::runtime_error("diagHermitian: d must be contiguous");
#endif

    //Set U = -M so eigenvalues will be sorted from largest to smallest
    if(isContiguous(M)) detail::copyNegElts(M.data(),makeRef(U));
    else                detail::copyNegElts(M.cbegin(),makeRef(U));


    auto info = detail::hermitianDiag(N,U.data(),d.data());
    if(info != 0) 
        {
        //println("M = \n",M);
        throw std::runtime_error("Error condition in diagHermitian");
        }

    //Correct the signs of the eigenvalues:
    d *= -1;
    //If M is transposed, we actually just computed the decomposition of
    //M^T=M^*, so conjugate U to compensate for this
    if(isTransposed(M)) conjugate(U);
    }

template<typename V>
void
diagGeneralRef(MatRefc<V> const& M,
               MatrixRef const& Rr,
               MatrixRef const& Ri,
               MatrixRef const& Lr,
               MatrixRef const& Li,
               VectorRef const& dr,
               VectorRef const& di);

template<class MatM, 
         class MatV,
         class Vecd,
         class>
void
eigen(MatM && M,
      MatV && Vr,
      MatV && Vi,
      Vecd && dr,
      Vecd && di)
    {
    auto N = ncols(M);
    resize(Vr,N,N);
    resize(Vi,N,N);
    resize(dr,N);
    resize(di,N);
    auto Lr = MatrixRef{};
    auto Li = MatrixRef{};
    diagGeneralRef(makeRef(M),makeRef(Vr),makeRef(Vi),Lr,Li,makeRef(dr),makeRef(di));
    }

template<class MatM, class MatV,class Vecd,class>
void
eigDecomp(MatM && M,
          MatV && Lr,
          MatV && Li,
          Vecd && dr,
          Vecd && di,
          MatV && Rr,
          MatV && Ri)
    {
    auto N = ncols(M);
    resize(Lr,N,N);
    resize(Li,N,N);
    resize(Rr,N,N);
    resize(Ri,N,N);
    resize(dr,N);
    resize(di,N);
    diagGeneralRef(makeRef(M),
                   makeRef(Rr),makeRef(Ri),
                   makeRef(Lr),makeRef(Li),
                   makeRef(dr),makeRef(di));
    }


template<typename T>
void
SVDRef(MatRefc<T> const& M,
       MatRef<T>  const& U, 
       VectorRef  const& D, 
       MatRef<T>  const& V,
       const Args & args);

template<class MatM, 
         class MatU,
         class VecD,
         class MatV,
         class>
void
SVD(MatM && M,
    MatU && U, 
    VecD && D, 
    MatV && V,
    const Args & args)
    {
    auto Mr = nrows(M),
         Mc = ncols(M);
    auto nsv = std::min(Mr,Mc);
    resize(U,Mr,nsv);
    resize(V,Mc,nsv);
    resize(D,nsv);
    SVDRef(makeRef(M),makeRef(U),makeRef(D),makeRef(V), args);
    }

template<class MatM, 
         class ScalarT,
         class>
Mat<common_type<val_type<MatM>,ScalarT>>
expHermitian(MatM && M,
             ScalarT t)
    {
    //using Mval = typename stdx::decay_t<MatM>::value_type;
    using valM = val_type<MatM>;
	Mat<valM> U;
    Vec<Real> d;
    diagHermitian(M,U,d);

    auto N = ncols(M);
	Mat<ScalarT> D(N,N);
    for(auto j : range(N))
        {
        D(j,j) = exp(d(j)*t);
        }
    auto expM = U*D*conj(transpose(U));
    return expM;
    }

namespace exptH_detail {

    int
    expPade(MatRef<Real> const& F, int N, int ideg);
    int
    expPade(MatRef<Cplx> const& F, int N, int ideg);

} //exptH_detail

template<class MatM,
         class ScalarT,
         class>
Mat<common_type<val_type<MatM>,ScalarT>>
expMatrix(MatM && M,
          ScalarT t,
          int ideg)
    {
    if(ideg <= 0) Error("PadeApproxDeg cannot be less than 1");
    auto N = ncols(M);
    if(N < 1) throw std::runtime_error("exp: 0 dimensional matrix");
    if(N != nrows(M))
        {
        printfln("M is %dx%d",nrows(M),ncols(M));
        throw std::runtime_error("exp: Input Matrix must be square");
        }
 
    auto tM = t*M;

    int info = 0;
    info = exptH_detail::expPade(makeRef(tM),N,ideg);

    if(info != 0)
        throw std::runtime_error("exp failed");

    return tM;
    }

  template<class MatA, 
         class MatQ,
         class MatR>
void
QR( MatA&& A,
    MatQ && Q,
    MatR && R,
    const Args & args)
   {

     auto complete = args.getBool("Complete", false);
     auto positive = args.getBool("PositiveDiagonal", false);
    using Aval = typename stdx::decay_t<MatA>::value_type;
    using Qval = typename stdx::decay_t<MatQ>::value_type;
    static_assert((isReal<Aval>() && isReal<Qval>()) || (isCplx<Aval>() && isCplx<Qval>()),
                  "A and Q must be both real or both complex in QR");
    int M = nrows(A);
    int N = ncols(A);
    if(M < 1 or N < 1) throw std::runtime_error("QR: 0 dimensional matrix");
    int Rrows = M;
    if (N > M)
      {
	println("Warning: QR for ncol > nrow may require pivoting to be numerically stable.");
	resize(Q,M,N);
      }
    else
      {
	Rrows = complete ? M : N;
	resize(Q,M,Rrows);
      }
    resize(R, Rrows, N);
  
#ifdef DEBUG
    if(!isContiguous(Q))
        throw std::runtime_error("QR: Q must be contiguous");
    if(!isContiguous(R))
        throw std::runtime_error("QR: R must be contiguous");
#endif
    for (int i = 0; i < M; ++i)
      {
      for (int j = 0; j < N; ++j)
	Q(i, j) = A(i,j);
      for (unsigned int j = N; j < ncols(Q); ++j)
	Q(i, j) = 0;
      }
    
    auto info = detail::QR(M,N,Rrows,Q.data(),R.data());
    if(info != 0) 
        {
        throw std::runtime_error("Error condition in QR");
        }
    if(N > M)
      reduceCols(Q,M);
    if (positive)
      {
	const int diagSize = Rrows < N ? Rrows : N;
	for (int i = 0; i < diagSize; i++)
	  {
	  if(std::real(R(i,i)) < 0)
	    {
	      row(R, i) *= -1.;
	      column(Q,i) *= -1.;
	    }
	 
	  }
      }
   }


  template<typename T>
  void
SVDRefLAPACK(
            MatRefc<T> const& M,
            MatRef<T>  const& U, 
            VectorRef  const& D, 
            MatRef<T>  const& V,
	     const Args & args)
{
    int Mr = nrows(M), 
         Mc = ncols(M);

    auto svdMethod = args.getString("SVDMethod", "automatic");

    auto pA = M.data();
    std::vector<T> cpA;
    cpA.resize(Mr*Mc);
    
    // LAPACK ?gesdd will read input matrix in column-major order. If we actually
    // want to perform SVD of M**T where M is stored in column-major, we have to pass
    // M**T stored in column-major. Copy of inpput matrix has to be done in any case, 
    // since input matrix is destroyed in ?gesdd
    if(isTransposed(M)) {
        for (unsigned int i=0; i<cpA.size(); i++, pA++) cpA[(i%Mc)*Mr + i/Mc] = *pA;
    } else {
        std::copy(pA,pA+Mr*Mc,cpA.data());
    }

    int info = -1;
    if (svdMethod == "automatic")
      {
      info = detail::SVD_gesdd(Mr, Mc, cpA.data(), U.data(), D.data(), V.data());

      // if gesdd failed, try gesvd; need to restore cpA data since gesdd destroyed it
      if(info != 0)
        {
          if(isTransposed(M)) {
              for (unsigned int i=0; i<cpA.size(); i++, pA++) cpA[(i%Mc)*Mr + i/Mc] = *pA;
          } else {
              std::copy(pA,pA+Mr*Mc,cpA.data());
          }
          info = detail::SVD_gesvd(Mr, Mc, cpA.data(), U.data(), D.data(), V.data());
        }
      }
    else if (svdMethod == "gesdd")
      info = detail::SVD_gesdd(Mr, Mc, cpA.data(), U.data(), D.data(), V.data());

    else if (svdMethod == "gesvd")
      info = detail::SVD_gesvd(Mr, Mc, cpA.data(), U.data(), D.data(), V.data());
    
    if(info != 0) 
      {
        throw std::runtime_error("Error condition in LAPACK SVD");
      }

    // from ?gesdd:
    // if JOBZ = 'S', V contains the first min(M=Mr,N=Mc) rows of
    // V**T (the right singular vectors, stored rowwise); 
    // Lapack stores V in column-major format, while the return of this function
    // expects row-major format of V, hence the V is reordered accordingly
    auto ncV = const_cast<T*>(V.data()); 
    auto pV  = reinterpret_cast<T*>(ncV);

    int l = std::min(Mr,Mc);
    std::vector<T> vt(l*Mc);
    std::copy(V.data(), V.data()+l*Mc, vt.data());
    for (unsigned int i=0; i<vt.size(); i++, pV++) *pV = detail::conjIfCplx(vt[(i%Mc)*l + i/Mc]);
    
}


} //namespace itensor

#endif
