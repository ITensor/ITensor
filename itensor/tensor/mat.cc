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
#include <limits>
#include "itensor/util/iterate.h"
#include "itensor/util/timers.h"
#include "itensor/tensor/lapack_wrap.h"
#include "itensor/tensor/mat.h"
#include "itensor/tensor/slicemat.h"
#include "itensor/detail/algs.h"

namespace itensor {

template<typename Iter1, typename Iter2, typename Func>
void
apply(Iter1 it1,
      Iter1 end1,
      Iter2 it2,
      Func const& f)
    {
    for(; it1 != end1; ++it1, ++it2)
        {
        f(*it1,*it2);
        }
    }

template<typename V, typename Func, typename Iter>
void
apply(MatRef<V> const& M,
      Iter it,
      Func const& f)
    {
    for(auto& el : M) 
        {
        f(el,*it);
        ++it;
        }
    }

void 
operator&=(MatrixRef const& a, MatrixRefc const& b)
    {
#ifdef DEBUG
    if(!(nrows(b)==nrows(a) && ncols(b)==ncols(a))) 
        throw std::runtime_error("mismatched sizes in MatrixRef operator&=");
#endif
    auto assign = [](Real& x, Real y) { x = y; };
    if(a.range()==b.range() && isContiguous(b))
        {
        auto pa = MAKE_SAFE_PTR(a.data(),a.store().size());
        auto pae = MAKE_SAFE_PTR_OFFSET(a.data(),dim(a.range()),a.store().size());
        auto pb = MAKE_SAFE_PTR(b.data(),b.store().size());
        apply(pa,pae,pb,assign);
        }
    else
        {
        apply(a,b.cbegin(),assign);
        }
    }

void
operator&=(CMatrixRef const& a, MatrixRefc const& b)
    {
#ifdef DEBUG
    if(!(nrows(b)==nrows(a) && ncols(b)==ncols(a)))
        throw std::runtime_error("mismatched sizes in MatrixRef operator&=");
#endif
    auto assign = [](Cplx& x, Real y) { x = y; };
    if(a.range()==b.range() && isContiguous(b))
        {
        auto pa = MAKE_SAFE_PTR(a.data(),a.store().size());
        auto pae = MAKE_SAFE_PTR_OFFSET(a.data(),dim(a.range()),a.store().size());
        auto pb = MAKE_SAFE_PTR(b.data(),b.store().size());
        apply(pa,pae,pb,assign);
        }
    else
        {
        apply(a,b.cbegin(),assign);
        }
    }

template<typename V>
void 
multReal(MatRef<V> const& M, Real fac)
    {
    if(isContiguous(M))
        {
#ifdef DEBUG
        if(M.size() > std::numeric_limits<unsigned long>::max()) 
            throw std::runtime_error("MatrixRef overflow of size beyond long unsigned int range");
#endif
        auto d = realData(M);
        dscal_wrapper(d.size(),fac,d.data());
        }
    else
        {
        for(auto& el : M) el *= fac;
        }
    }

void
multCplx(CMatrixRef const& M, Cplx fac)
    {
    for(auto& el : M) el *= fac;
    } 

void 
operator*=(MatrixRef const& M, Real fac)
    {
    multReal(M,fac);
    }
void 
operator*=(CMatrixRef const& M, Real fac)
    {
    multReal(M,fac);
    }

void
operator*=(CMatrixRef const& M, Cplx fac)
    {
    multCplx(M,fac);
    }

template<typename V>
void 
divReal(MatRef<V> const& M, Real fac)
    {
    if(isContiguous(M))
        {
        auto d = realData(M);
        auto mend = d.data()+d.size();
        for(auto m = d.data(); m != mend; ++m)
            {
            *m /= fac;
            }
        }
    else
        {
        for(auto& el : M) el /= fac;
        }
    }

void 
operator/=(MatrixRef const& M, Real fac)
    {
    divReal(M,fac);
    }
void 
operator/=(CMatrixRef const& M, Real fac)
    {
    divReal(M,fac);
    }

template<typename MatT1, typename MatT2>
void
call_daxpy(MatT1& A, MatT2 const& B, Real alpha_)
    {
    LAPACK_REAL alpha = alpha_;
    LAPACK_INT inc = 1;
    auto Ad = realData(A);
    auto Bd = realData(B);
#ifdef DEBUG
    if(Ad.size() != Bd.size())
        throw std::runtime_error("mismatched sizes in MatrixRef/Matrix call_daxpy");
    if(Ad.size() > std::numeric_limits<unsigned long>::max()) 
        throw std::runtime_error("overflow of size beyond long unsigned int range");
#endif
    daxpy_wrapper(Ad.size(),alpha,Bd.data(),inc,Ad.data(),inc);
    }

template<typename V>
void
add(MatRef<V> const& A, MatRefc<V> const& B)
    {
#ifdef DEBUG
    if(!(ncols(A)==ncols(B) && nrows(A)==nrows(B))) 
        throw std::runtime_error("MatrixRef +=: mismatched sizes");
#endif
    if(B.range()==A.range() && isContiguous(B))
        {
        call_daxpy(A,B,+1);
        }
    else
        {
        auto pluseq = [](V& x, V y) { x += y; };
        apply(A,B.cbegin(),pluseq);
        }
    }
void
operator+=(MatrixRef const& A, MatrixRefc const& B)
    {
    add(A,B);
    }
void
operator+=(MatrixRef const& A, Matrix && B)
    {
    add(A,makeRef(B));
    }
void
operator+=(CMatrixRef const& A, CMatrixRefc const& B)
    {
    add(A,B);
    }
void
operator+=(CMatrixRef const& A, CMatrix && B)
    {
    add(A,makeRef(B));
    }


template<typename V>
void
subtract(MatRef<V> const& A, MatRefc<V> const& B)
    {
#ifdef DEBUG
    if(!(ncols(A)==ncols(B) && nrows(A)==nrows(B))) 
        throw std::runtime_error("MatrixRef +=: mismatched sizes");
#endif
    if(B.range()==A.range() && isContiguous(B))
        {
        call_daxpy(A,B,-1);
        }
    else
        {
        auto minuseq = [](V& x, V y) { x -= y; };
        apply(A,B.cbegin(),minuseq);
        }
    }
void
operator-=(MatrixRef const& A, MatrixRefc const& B)
    {
    subtract(A,B);
    }
void
operator-=(MatrixRef const& A, Matrix && B)
    {
    subtract(A,makeRef(B));
    }
void
operator-=(CMatrixRef const& A, CMatrixRefc const& B)
    {
    subtract(A,B);
    }
void
operator-=(CMatrixRef const& A, CMatrix && B)
    {
    subtract(A,makeRef(B));
    }

template<typename V>
void
printMatrix(std::ostream& s, MatRefc<V> const& M)
    {
    for(auto r : range(nrows(M)))
        {
        s << "|";
        for(auto c : range(ncols(M)))
            {
            s << formatVal(M(r,c));
            s << (1+c == ncols(M) ? "|" : " ");
            }
        if(r < nrows(M)) s << "\n";
        }
    }

template<>
std::ostream&
operator<<(std::ostream& s, MatrixRefc const& M)
    {
    printMatrix(s,M);
    return s;
    }
template<>
std::ostream&
operator<<(std::ostream& s, CMatrixRefc const& M)
    {
    printMatrix(s,M);
    return s;
    }


//#ifdef METHOD2
//    if(not (alpha == Cplx(1.,0.) && beta == Cplx(0.,0.)))
//        {
//        throw std::runtime_error("alpha,beta must be 1,0");
//        }
//
//    auto aCsize = m*k;
//    auto bCsize = k*n;
//    auto cCsize = m*n;
//    //WARNING: logically const, not thread safe!
//    auto Ad = reinterpret_cast<Real*>(const_cast<Cplx*>(A));
//    auto Bd = reinterpret_cast<Real*>(const_cast<Cplx*>(B));
//    auto Cd = reinterpret_cast<Real*>(C);
//    auto arb = Ad;
//    auto aib = Ad+aCsize;
//    auto brb = Bd;
//    auto bib = Bd+bCsize;
//    auto crb = Cd;
//    auto cib = Cd+cCsize;
//        //print("re(A):");
//        //for(auto i = arb; i != aib; ++i)
//        //    {
//        //    print(" ",*i);
//        //    }
//        //println();
//        //print("im(A):");
//        //for(auto i = aib; i != aib+aCsize; ++i)
//        //    {
//        //    print(" ",*i);
//        //    }
//        //println();
//    toCplx(Ad,aCsize,false);
//    toCplx(Bd,bCsize,false);
//
//    //print("re(A):");
//    //for(auto i = arb; i != aib; ++i)
//    //    {
//    //    print(" ",*i);
//    //    }
//    //println();
//    //print("im(A):");
//    //for(auto i = aib; i != aib+aCsize; ++i)
//    //    {
//    //    print(" ",*i);
//    //    }
//    //println();
//
//    cblas_dgemm(CblasColMajor,at,bt,m,n,k,+1,arb,lda,brb,ldb,0.,crb,m);
//    cblas_dgemm(CblasColMajor,at,bt,m,n,k,-1,aib,lda,bib,ldb,1.,crb,m);
//    cblas_dgemm(CblasColMajor,at,bt,m,n,k,+1,arb,lda,bib,ldb,0.,cib,m);
//    cblas_dgemm(CblasColMajor,at,bt,m,n,k,+1,aib,lda,brb,ldb,1.,cib,m);
//    toCplx(Ad,aCsize,true);
//    toCplx(Bd,bCsize,true);
//    toCplx(Cd,cCsize,true);
//#endif //METHOD2








//void
//mult(CMatrixRefc A,
//     CMatrixRefc B,
//     CMatrixRef  C)
//    {
//    gemm(A,B,C,1.,0.);
//    }

void
multAdd(MatrixRefc A, 
        MatrixRefc B, 
        MatrixRef  C)
    {
    gemm(A,B,C,1.,1.);
    }


template<typename V>
void
call_gemv(MatRefc<V> const& M,
          VecRefc<V> const& x, 
          VecRef<V>       & y,
          Real alpha,
          Real beta,
          bool fromleft)
    {
#ifdef DEBUG
    if(!isContiguous(M))
        throw std::runtime_error("multiplication of non-contiguous matrixref by vector not currently supported");
#endif
    auto trans = fromleft;
    LAPACK_INT m = nrows(M),
               n = ncols(M);
    if(isTransposed(M))
        {
        trans = !fromleft;
        m = ncols(M);
        n = nrows(M);
        }
    gemv_wrapper(trans,alpha,beta,m,n,M.data(),x.data(),stride(x),y.data(),stride(y));
    }
template void call_gemv(MatRefc<Real> const&,VecRefc<Real> const&, VecRef<Real>&,Real,Real,bool);
template void call_gemv(MatRefc<Cplx> const&,VecRefc<Cplx> const&, VecRef<Cplx>&,Real,Real,bool);


template<typename V>
void
mult(MatRefc<V> M,
     VecRefc<V> x,
     VecRef<V> y,
     bool fromleft)
    {
#ifdef DEBUG
    if(fromleft ? nrows(M)!=x.size() : ncols(M)!=x.size()) 
        throw std::runtime_error("matrix vector mult: mismatched sizes");
    if(fromleft ? ncols(M)!=y.size() : nrows(M)!=y.size())
        throw std::runtime_error("matrix vector mult: wrong size for result (y) vec");
#endif
    call_gemv(M,x,y,1,0,fromleft);
    }
template void mult(MatRefc<Real>,VecRefc<Real>,VecRef<Real>,bool);
template void mult(MatRefc<Cplx>,VecRefc<Cplx>,VecRef<Cplx>,bool);

template<typename VM, typename Vx>
Vec<common_type<VM,Vx>>
mult(MatRefc<VM> M,
     VecRefc<Vx> x,
     bool fromleft)
    {
    auto res = Vec<common_type<VM,Vx>>(x.size());
    mult(M,x,makeRef(res));
    return res;
    }
template Vector mult(MatrixRefc,VectorRefc, bool);
//template CVector mult(CMatrixRefc,VectorRefc, bool);
//template CVector mult(MatrixRefc,CVectorRefc, bool);
template CVector mult(CMatrixRefc,CVectorRefc, bool);

template<typename V>
void
multAdd(MatRefc<V> M,
        VecRefc<V> x,
        VecRef<V> y,
        bool fromleft)
    {
#ifdef DEBUG
    if(fromleft ? nrows(M)!=x.size() : ncols(M)!=x.size()) 
        throw std::runtime_error("multAdd: mismatched sizes");
    if(fromleft ? ncols(M)!=y.size() : nrows(M)!=y.size())
        throw std::runtime_error("multAdd: wrong size for result (y) vec");
#endif
    call_gemv(M,x,y,1,1,fromleft);
    }
template void multAdd(MatRefc<Real>,VecRefc<Real>,VecRef<Real>,bool);
template void multAdd(MatRefc<Cplx>,VecRefc<Cplx>,VecRef<Cplx>,bool);


template<typename V>
void
multSub(MatRefc<V> M,
        VecRefc<V> x,
        VecRef<V> y,
        bool fromleft)
    {
#ifdef DEBUG
    if(fromleft ? nrows(M)!=x.size() : ncols(M)!=x.size()) 
        throw std::runtime_error("multSub: mismatched sizes");
    if(fromleft ? ncols(M)!=y.size() : nrows(M)!=y.size())
        throw std::runtime_error("multSub: wrong size for result (y) vec");
#endif
    call_gemv(M,x,y,-1,1,fromleft);
    }
template void multSub(MatRefc<Real>,VecRefc<Real>,VecRef<Real>,bool);
template void multSub(MatRefc<Cplx>,VecRefc<Cplx>,VecRef<Cplx>,bool);

Matrix
eye(size_t Nr, size_t Nc)
    {
    auto M = Matrix(Nr,Nc);
    auto N = std::min(Nr,Nc);
    for(auto j : range(N)) M(j,j) = 1.0;
    return M;
    }

} //namespace itensor
