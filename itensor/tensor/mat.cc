//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#include <limits>
#include "itensor/util/count.h"
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
        auto pae = MAKE_SAFE_PTR_OFFSET(a.data(),area(a.range()),a.store().size());
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
        if(M.size() > std::numeric_limits<LAPACK_INT>::max()) 
            throw std::runtime_error("MatrixRef overflow of size beyond LAPACK_INT range");
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
operator*=(MatrixRef const& M, Real fac)
    {
    multReal(M,fac);
    }
void 
operator*=(CMatrixRef const& M, Real fac)
    {
    multReal(M,fac);
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
    if(Ad.size() > std::numeric_limits<LAPACK_INT>::max()) 
        throw std::runtime_error("overflow of size beyond LAPACK_INT range");
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

void
randomize(MatrixRef const& M)
    {
    for(auto& el : M) el = detail::quickran();
    }
void
randomize(Matrix & M)
    {
    for(auto& el : M) el = detail::quickran();
    }
void
randomize(CMatrixRef const& M)
    {
    for(auto& el : M) el = Cplx(detail::quickran(),detail::quickran());
    }
void
randomize(CMatrix & M)
    {
    for(auto& el : M) el = Cplx(detail::quickran(),detail::quickran());
    }


template<typename V>
void
printMatrix(std::ostream& s, MatRefc<V> const& M)
    {
    for(auto r : count(nrows(M)))
        {
        s << "|";
        for(auto c : count(ncols(M)))
            {
            s << detail::printVal(M(r,c));
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

// C = alpha*A*B + beta*C
template<typename V>
void
gemm(MatRefc<V> A, 
     MatRefc<V> B, 
     MatRef<V>  C,
     Real alpha,
     Real beta)
    {
#ifdef DEBUG
    if(!(isContiguous(A) && isContiguous(B) && isContiguous(C))) 
        throw std::runtime_error("multiplication of non-contiguous MatrixRefs not currently supported");
#endif
    if(isTransposed(C))
        {
        //Do C = Bt*At instead of Ct=A*B
        //Recall that C.data() points to elements of C, not C.t()
        //regardless of whether C.transpose()==true or false
        std::swap(A,B);
        A = transpose(A);
        B = transpose(B);
        C = transpose(C);
        }

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
    START_TIMER(33)
    gemm_wrapper(isTransposed(A),isTransposed(B),
                 nrows(A),ncols(B),ncols(A),
                 alpha,A.data(),B.data(),beta,C.data());
    STOP_TIMER(33)
    }

void
mult(MatrixRefc A, 
     MatrixRefc B, 
     MatrixRef  C)
    {
    gemm(A,B,C,1.,0.);
    }

void
multAdd(MatrixRefc A, 
        MatrixRefc B, 
        MatrixRef  C)
    {
    gemm(A,B,C,1.,1.);
    }

void
mult(CMatrixRefc A,
     CMatrixRefc B,
     CMatrixRef  C)
    {
    gemm(A,B,C,1.,0.);
    }

void
call_dgemv(MatrixRefc const& M,
           VectorRefc const& x, 
           VectorRef       & y,
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
    dgemv_wrapper(trans,alpha,beta,m,n,M.data(),x.data(),stride(x),y.data(),stride(y));
    }

void
mult(MatrixRefc M,
     VectorRefc x,
     VectorRef y,
     bool fromleft)
    {
#ifdef DEBUG
    if(fromleft ? nrows(M)!=x.size() : ncols(M)!=x.size()) 
        throw std::runtime_error("matrix vector mult: mismatched sizes");
    if(fromleft ? ncols(M)!=y.size() : nrows(M)!=y.size())
        throw std::runtime_error("matrix vector mult: wrong size for result (y) vec");
#endif
    call_dgemv(M,x,y,1,0,fromleft);
    }

void
multAdd(MatrixRefc M,
        VectorRefc x,
        VectorRef y,
        bool fromleft)
    {
#ifdef DEBUG
    if(fromleft ? nrows(M)!=x.size() : ncols(M)!=x.size()) 
        throw std::runtime_error("multAdd: mismatched sizes");
    if(fromleft ? ncols(M)!=y.size() : nrows(M)!=y.size())
        throw std::runtime_error("multAdd: wrong size for result (y) vec");
#endif
    call_dgemv(M,x,y,1,1,fromleft);
    }

void
multSub(MatrixRefc M,
        VectorRefc x,
        VectorRef y,
        bool fromleft)
    {
#ifdef DEBUG
    if(fromleft ? nrows(M)!=x.size() : ncols(M)!=x.size()) 
        throw std::runtime_error("multSub: mismatched sizes");
    if(fromleft ? ncols(M)!=y.size() : nrows(M)!=y.size())
        throw std::runtime_error("multSub: wrong size for result (y) vec");
#endif
    call_dgemv(M,x,y,-1,1,fromleft);
    }


} //namespace itensor
