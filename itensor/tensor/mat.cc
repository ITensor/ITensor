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

template<typename Func, typename Iter>
void
apply(MatrixRef const& v,
      Iter it,
      Func const& f)
    {
    for(auto& el : v) 
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

void 
operator*=(MatrixRef const& A, Real fac)
    {
    if(isContiguous(A))
        {
#ifdef DEBUG
        if(A.size() > std::numeric_limits<LAPACK_INT>::max()) 
            throw std::runtime_error("MatrixRef overflow of size beyond LAPACK_INT range");
#endif
        dscal_wrapper(A.size(),fac,A.data());
        }
    else
        {
        for(auto& el : A) el *= fac;
        }
    }

void 
operator/=(MatrixRef const& A, Real fac)
    {
    //if(fac == 0) throw std::runtime_error("MatrixRef /=: divide by zero");
    //operator*=(A,1./fac);
    if(isContiguous(A))
        {
        auto ae = A.data()+A.size();
        for(auto a = A.data(); a != ae; ++a)
            {
            *a /= fac;
            }
        }
    else
        {
        for(auto& el : A) el /= fac;
        }
    }

template<typename MatT1, typename MatT2>
void
call_daxpy(MatT1& A, MatT2 const& B, Real alpha_)
    {
    LAPACK_REAL alpha = alpha_;
    LAPACK_INT inc = 1;
    LAPACK_INT size = A.size();
#ifdef DEBUG
    if(A.size() != B.size())
        throw std::runtime_error("mismatched sizes in MatrixRef/Matrix call_daxpy");
    if(A.size() > std::numeric_limits<LAPACK_INT>::max()) 
        throw std::runtime_error("overflow of size beyond LAPACK_INT range");
#endif
    daxpy_wrapper(size,alpha,B.data(),inc,A.data(),inc);
    }

void
operator+=(MatrixRef const& A, MatrixRefc const& B)
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
        auto pluseq = [](Real& x, Real y) { x += y; };
        apply(A,B.cbegin(),pluseq);
        }
    }

void
operator-=(MatrixRef const& A, MatrixRefc const& B)
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
        auto minuseq = [](Real& x, Real y) { x -= y; };
        apply(A,B.cbegin(),minuseq);
        }
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


template<>
std::ostream&
operator<<(std::ostream& s, MatrixRefc const& M)
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
    return s;
    }

template<>
std::ostream&
operator<<(std::ostream& s, CMatrixRefc const& M)
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
    return s;
    }

// C = alpha*A*B + beta*C
template<typename V>
void
call_gemm(MatRefc<V> A, 
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
    gemm_wrapper(isTransposed(A),isTransposed(B),
                 nrows(A),ncols(B),ncols(A),
                 alpha,A.data(),B.data(),beta,C.data());
    }

void
mult(MatrixRefc A, 
     MatrixRefc B, 
     MatrixRef  C)
    {
    call_gemm(A,B,C,1.,0.);
    }

void
multAdd(MatrixRefc A, 
        MatrixRefc B, 
        MatrixRef  C)
    {
    call_gemm(A,B,C,1.,1.);
    }

void
mult(CMatrixRefc A,
     CMatrixRefc B,
     CMatrixRef  C)
    {
    call_gemm(A,B,C,1.,0.);
    }

void
reduceCols(Matrix & M, size_t new_ncols)
    {
#ifdef DEBUG
    if(new_ncols > ncols(M)) throw std::runtime_error("new ncols > old ncols in reduceCols");
#endif
    M.resize(MatRange(nrows(M),new_ncols));
    }

void
resize(Matrix & M, size_t nrows, size_t ncols)
    {
    M.resize(MatRange(nrows,ncols));
    }

void
call_dgemv(const MatrixRefc& M,
           const VectorRefc& x, 
           VectorRef& y,
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

Matrix 
operator*(MatrixRefc const& A, Real fac)
    { 
    Matrix res(A);
    res *= fac; 
    return res; 
    }

Matrix 
operator*(Real fac, MatrixRefc const& A)
    { 
    Matrix res(A);
    res *= fac; 
    return res; 
    }

Matrix 
operator*(Matrix && A, Real fac)
    { 
    Matrix res(std::move(A));
    res *= fac; 
    return res; 
    }

Matrix 
operator*(Real fac, Matrix && A)
    { 
    Matrix res(std::move(A));
    res *= fac; 
    return res; 
    }

Matrix 
operator/(MatrixRefc const& A, Real fac)
    { 
    Matrix res(A);
    res /= fac; 
    return res; 
    }

Matrix 
operator/(Matrix && A, Real fac)
    { 
    Matrix res(std::move(A));
    res /= fac; 
    return res; 
    }

Matrix 
operator+(MatrixRefc const& A, MatrixRefc const& B)
    { 
    Matrix res(A);
    res += B; 
    return res; 
    }

Matrix 
operator+(MatrixRefc const& A, Matrix && B) 
    { 
    Matrix res(std::move(B)); 
    res += A; 
    return res; 
    }

Matrix 
operator+(Matrix && A, MatrixRefc const& B) 
    { 
    Matrix res(std::move(A)); 
    res += B; 
    return res; 
    }

Matrix 
operator+(Matrix && A, Matrix && B)
    {
    Matrix mA(std::move(A)); 
    Matrix mB(std::move(B)); 
    mA += mB; 
    return mA; 
    }

Matrix 
operator-(MatrixRefc const& A, MatrixRefc const& B)
    { 
    Matrix res(A);
    res -= B; 
    return res; 
    }

Matrix 
operator-(MatrixRefc const& A, Matrix && B) 
    { 
    Matrix res(std::move(B)); 
    res *= -1;
    res += A; 
    return res; 
    }

Matrix 
operator-(Matrix && A, MatrixRefc const& B) 
    { 
    Matrix res(std::move(A)); 
    res -= B; 
    return res; 
    }

Matrix 
operator-(Matrix && A, Matrix && B)
    {
    Matrix mA(std::move(A)); 
    Matrix mB(std::move(B)); 
    mA -= mB; 
    return mA; 
    }

Matrix 
matrixMult(MatrixRefc const& A,
           MatrixRefc const& B)
    {
    Matrix C(nrows(A),ncols(B));
    call_gemm(A,B,makeRef(C),1.,0.);
    return C;
    }

Vector
operator*(MatrixRefc const& A,
          VectorRefc const& v)
    {
    Vector res(nrows(A));
    mult(A,v,res);
    return res;
    }

Vector
operator*(VectorRefc const& v,
          MatrixRefc const& A)
    {
    Vector res(ncols(A));
    bool fromleft = true;
    mult(A,v,res,fromleft);
    return res;
    }

} //namespace itensor
