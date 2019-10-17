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

namespace exptH_detail {
    //Helper for expPade
    template<typename Iter, typename V>
    void
    multTElts(Iter m,
              MatRef<V> const& F,
              Cplx t)	
        {
        auto fe = F.data()+F.size();
        if(t.imag() == 0)
            {
            Real tre = t.real();
            for(auto f = F.data(); f != fe; ++f, ++m)
                {
                *f = tre * (*m);
                }
            }
        else
            {
            for(auto f = F.data(); f != fe; ++f, ++m) 
                {
                *f = t * (*m);
                }
            }
        }
	
    //Helper for expPade
    template<typename Iter>
    void
    multTElts(Iter mre,
              Iter mim,
              std::vector<Cplx> & Hc,
              Cplx t)
        {
        if(t.imag() == 0)
            {
            for(auto& z : Hc)
                {
                realRef(z) = t.real() * (*mre);
                imagRef(z) = t.real() * (*mim);
                ++mre;
                ++mim;
                }
            }
        else
            {
            for(auto& z : Hc)
                {
                realRef(z) = t.real() * (*mre) - t.imag() * (*mim);
                imagRef(z) = t.real() * (*mim) + t.imag() * (*mre);
                ++mre;
                ++mim;
                }
            }
        }
	
    //template <typename T>
    //int padeExp(VecRef<T> const& y, int N, MatRef<T> const& F, int ideg);
    int
    padeExp(MatRef<Real> const& F, int N, int ideg);
    int
    padeExp(MatRef<Cplx> const& F, int N, int ideg);

} //exptH_detail

template<class Vecy,
         class MatH,
         class>
void
expMatrixApply(Vecy && y,
               MatH && H,
               Cplx t,
               int ideg)
    {
    if(ideg <= 0)
        {
        Error("ideg cannot be less than 1");
        }
    else
        {
        using yval = typename stdx::decay_t<Vecy>::value_type;
        using Hval = typename stdx::decay_t<MatH>::value_type;
        if(t.imag() == 0)
            static_assert( !(isReal<yval>() && isCplx<Hval>()),
                          "expMatrixApply: real t - y must complex when H is complex");	
        else
            static_assert( isCplx<yval>(),
                          "expMatrixApply: cplx t - y must be complex");

        auto N = ncols(H);
        if(N < 1) throw std::runtime_error("expMatrixApply: 0 dimensional matrix");
        if(N != nrows(H))
            {
            printfln("H is %dx%d",nrows(H),ncols(H));
            throw std::runtime_error("expMatrixApply: Input Matrix must be square");
            }

        auto F = Mat<yval>(nrows(H),ncols(H));

#ifdef DEBUG
        if(!isContiguous(y))
            throw std::runtime_error("expMatrixApply: y must be contiguous");
#endif
        
        //TODO: use mult for matrix directly
        if(isContiguous(H)) 
            {
            exptH_detail::multTElts(H.data(),makeRef(F),t);
            }
        else                
            {
            exptH_detail::multTElts(H.cbegin(),makeRef(F),t);
            }

        int info = 0;
        info = exptH_detail::padeExp(makeRef(F),N,ideg);

        if(info != 0)
            {
            throw std::runtime_error("expMatrixApply failed");
            }

        y &= mult(makeRef(F),y);
        }
    }

} //namespace itensor

#endif
