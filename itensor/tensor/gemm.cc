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
#include "itensor/tensor/lapack_wrap.h"
#include "itensor/tensor/slicemat.h"
#include "itensor/util/safe_ptr.h"

namespace itensor {

struct dgemmTask
    {
    bool copyToC = false;
    bool copyFromC = false;
    size_t Apart = 0ul,
           Bpart = 0ul,
           Cpart = 0ul;
    Real alpha = 0.,
         beta  = 0.;

    dgemmTask(size_t Ap,
              size_t Bp,
              size_t Cp,
              Real a,
              Real b)
      : Apart(Ap),Bpart(Bp),Cpart(Cp),alpha(a),beta(b)
        { 
        if(b != 0.0) copyFromC = true;
        }

    dgemmTask(size_t Ap,
              size_t Bp,
              size_t Cp,
              Real a)
      : Apart(Ap),Bpart(Bp),Cpart(Cp),alpha(a),beta(1.0)
        { }

    dgemmTask(size_t Cp)
       : copyToC(true),Cpart(Cp)
        { }
    };

void
cplxToRealBuf(SAFE_PTR_OF(const Real) C,
              size_t imagPart,
              SAFE_PTR_OF(Real) Rbegin,
              SAFE_PTR_OF(Real) Rend)
    {
    C += imagPart;
    for(auto r=Rbegin; r != Rend; ++r, C+=2)
        {
        *r = *C;
        }
    }
void
realBufToCplx(SAFE_PTR_OF(Real) Rbegin,
              SAFE_PTR_OF(Real) Rend,
              SAFE_PTR_OF(Real) C,
              size_t imagPart)
    {
    C += imagPart;
    for(auto r=Rbegin; r != Rend; ++r, C+=2)
        {
        *C = *r;
        }
    }

template<size_t NTask, typename VA, typename VB>
void
gemm_emulator(MatRefc<VA> A,
              MatRefc<VB> B,
              MatRef<Cplx>  C,
              Real alpha,
              Real beta,
              std::array<const dgemmTask,NTask> const& tasks)
    {
    static_assert(!(isReal(A) && isReal(B)),
                  "only use gemm_emulator for complex case");
    auto Abufsize = isCplx(A) ? A.size() : 0ul;
    auto Bbufsize = isCplx(B) ? B.size() : 0ul;
    auto Cbufsize = C.size();

    auto Ad = MAKE_SAFE_PTR(A.data(),A.size());
    auto Bd = MAKE_SAFE_PTR(B.data(),B.size());
    auto Cd = MAKE_SAFE_PTR(C.data(),C.size());
    auto Ard = SAFE_REINTERPRET(const Real,Ad);
    auto Brd = SAFE_REINTERPRET(const Real,Bd);
    auto Crd = SAFE_REINTERPRET(Real,Cd);

    auto d = vector_no_init<Real>(Abufsize+Bbufsize+Cbufsize);
    auto pd = MAKE_SAFE_PTR(d.data(),d.size());
    auto ab = pd;
    auto ae = ab+Abufsize;
    auto bb = ae;
    auto be = bb+Bbufsize;
    auto cb = be;
    auto ce = cb+Cbufsize;

    for(auto t : tasks)
        {
        if(t.copyToC)
            {
            realBufToCplx(cb,ce,Crd,t.Cpart);
            }
        else
            {
            if(isCplx(A))   cplxToRealBuf(Ard,t.Apart,ab,ae);
            if(isCplx(B))   cplxToRealBuf(Brd,t.Bpart,bb,be);
            if(t.copyFromC) cplxToRealBuf(Crd,t.Cpart,cb,ce);

            if(isReal(A))
                {
                gemm_wrapper(isTransposed(A),
                             isTransposed(B),
                             nrows(A),
                             ncols(B),
                             ncols(A),
                             t.alpha,
                             SAFE_PTR_GET(Ard,A.size()),
                             SAFE_PTR_GET(bb,B.size()),
                             t.beta,
                             SAFE_PTR_GET(cb,C.size()));
                }
            else if(isReal(B))
                {
                gemm_wrapper(isTransposed(A),
                             isTransposed(B),
                             nrows(A),
                             ncols(B),
                             ncols(A),
                             t.alpha,
                             SAFE_PTR_GET(ab,A.size()),
                             SAFE_PTR_GET(Brd,B.size()),
                             t.beta,
                             SAFE_PTR_GET(cb,C.size()));
                }
            else //A and B are both Cplx
                {
                gemm_wrapper(isTransposed(A),
                             isTransposed(B),
                             nrows(A),
                             ncols(B),
                             ncols(A),
                             t.alpha,
                             SAFE_PTR_GET(ab,A.size()),
                             SAFE_PTR_GET(bb,B.size()),
                             t.beta,
                             SAFE_PTR_GET(cb,C.size()));
                }
            }
        }
    }

void
gemm_impl(MatRefc<Cplx> A,
          MatRefc<Cplx> B,
          MatRef<Cplx>  C,
          Real alpha,
          Real beta)
    {
#ifdef ITENSOR_USE_ZGEMM
    gemm_wrapper(isTransposed(A),
                 isTransposed(B),
                 nrows(A),
                 ncols(B),
                 ncols(A),
                 alpha,
                 A.data(),
                 B.data(),
                 beta,
                 C.data());
#else //emulate zgemm by calling dgemm four times
    std::array<const dgemmTask,6> 
    tasks = 
        {{dgemmTask(0,0,0,+alpha,beta),
          dgemmTask(1,1,0,-alpha),
          dgemmTask(0),
          dgemmTask(1,0,1,+alpha,beta),
          dgemmTask(0,1,1,+alpha),
          dgemmTask(1)
          }};
    gemm_emulator(A,B,C,alpha,beta,tasks);
#endif
    }


void
gemm_impl(MatRefc<Real> A,
          MatRefc<Cplx> B,
          MatRef<Cplx>  C,
          Real alpha,
          Real beta)
    {
    std::array<const dgemmTask,4> 
    tasks = 
        {{dgemmTask(0,0,0,+alpha,beta),
          dgemmTask(0),
          dgemmTask(0,1,1,+alpha,beta),
          dgemmTask(1)
          }};
    gemm_emulator(A,B,C,alpha,beta,tasks);
    }

void
gemm_impl(MatRefc<Cplx> A,
          MatRefc<Real> B,
          MatRef<Cplx>  C,
          Real alpha,
          Real beta)
    {
    std::array<const dgemmTask,4> 
    tasks = 
        {{dgemmTask(0,0,0,+alpha,beta),
          dgemmTask(0),
          dgemmTask(1,0,1,+alpha,beta),
          dgemmTask(1)
          }};
    gemm_emulator(A,B,C,alpha,beta,tasks);
    }

void
gemm_impl(MatRefc<Real> A,
          MatRefc<Real> B,
          MatRef<Real>  C,
          Real alpha,
          Real beta)
    {
    //call dgemm directly
    gemm_wrapper(isTransposed(A),
                 isTransposed(B),
                 nrows(A),
                 ncols(B),
                 ncols(A),
                 alpha,
                 A.data(),
                 B.data(),
                 beta,
                 C.data());
    }

// C = alpha*A*B + beta*C
template<typename VA, typename VB>
void
gemm(MatRefc<VA> A, 
     MatRefc<VB> B, 
     MatRef<common_type<VA,VB>>  C,
     Real alpha,
     Real beta)
    {
#ifdef DEBUG
    if(!(isContiguous(A) && isContiguous(B) && isContiguous(C))) 
        throw std::runtime_error("multiplication of non-contiguous MatrixRefs not currently supported");
#endif

#ifdef DEBUG
    if(ncols(A) != nrows(B))
        throw std::runtime_error("matrices A, B incompatible");
    if(nrows(A) != nrows(C) || ncols(B) != ncols(C))
        {
        printfln("A is %dx%d",nrows(A),ncols(A));
        printfln("B is %dx%d",nrows(B),ncols(B));
        printfln("C is %dx%d",nrows(C),ncols(C));
        throw std::runtime_error("mult(_add) AxB -> C: matrix C incompatible");
        }
#endif
    if(isTransposed(C))
        {
        //Do C = Bt*At instead of Ct=A*B
        //Recall that C.data() points to elements of C, not C.t()
        //regardless of whether C.transpose()==true or false
        gemm_impl(transpose(B),transpose(A),transpose(C),alpha,beta);
        }
    else
        {
        gemm_impl(A,B,C,alpha,beta);
        }
    }
template void gemm(MatRefc<Real>, MatRefc<Real>, MatRef<Real>,Real,Real);
template void gemm(MatRefc<Real>, MatRefc<Cplx>, MatRef<Cplx>,Real,Real);
template void gemm(MatRefc<Cplx>, MatRefc<Real>, MatRef<Cplx>,Real,Real);
template void gemm(MatRefc<Cplx>, MatRefc<Cplx>, MatRef<Cplx>,Real,Real);


} //namespace itensor
