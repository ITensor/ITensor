//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
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

    int
    hermitianDiag(int N, Real *Udata, Real *ddata);
    int
    hermitianDiag(int N, Cplx *Udata,Real *ddata);
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
       Real thresh);

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
    Real thresh)
    {
    auto Mr = nrows(M),
         Mc = ncols(M);
    auto nsv = std::min(Mr,Mc);
    resize(U,Mr,nsv);
    resize(V,Mc,nsv);
    resize(D,nsv);
    SVDRef(makeRef(M),makeRef(U),makeRef(D),makeRef(V),thresh);
    }

} //namespace itensor

#endif
