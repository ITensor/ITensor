//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#include "mat.h"
#include "lapack_wrap.h"
#include <limits>
#include "detail/algs.h"

namespace itensor {

//
//
//  .      .    __    _______   ____    .____  .____
//  |\    /|   /  \      |     |    |   |      |    
//  | \  / |  |____|     |     |___/    |___   |--- 
//  |  \/  |  |    |     |     |   \    |      |    
//  |      |  |    |     |     |    |   |____  |    
//
//

template<typename Func, typename Iter>
void
apply(const MatRef& v,
      Iter it,
      const Func& f)
    {
    for(auto& el : v) 
        {
        f(el,*it);
        ++it;
        }
    }

void 
operator&=(const MatRef& a, MatRefc b)
    {
#ifdef DEBUG
    if(!(b.Nrows()==a.Nrows() && b.Ncols()==a.Ncols())) 
        throw std::runtime_error("mismatched sizes in MatRef operator&=");
#endif
    auto assign = [](Real& x, Real y) { x = y; };
    if(a.ind()==b.ind() && b.contiguous())
        {
        apply(a,b.data(),assign);
        }
    else
        {
        apply(a,b.cbegin(),assign);
        }
    }

void 
operator*=(const MatRef& a, Real fac)
    {
    if(a.contiguous())
        {
#ifdef DEBUG
        if(a.size() > std::numeric_limits<LAPACK_INT>::max()) 
            throw std::runtime_error("MatRef overflow of size beyond LAPACK_INT range");
#endif
        dscal_wrapper(a.size(),fac,a.data());
        }
    else
        {
        for(auto& el : a) el *= fac;
        }
    }

void 
operator/=(const MatRef& a, Real fac)
    {
    if(fac == 0) throw std::runtime_error("MatRef /=: divide by zero");
    operator*=(a,1./fac);
    }

template<typename MatT1, typename MatT2>
void
call_daxpy(MatT1& A, const MatT2& B, Real alpha_)
    {
    LAPACK_REAL alpha = alpha_;
    LAPACK_INT inc = 1;
    LAPACK_INT size = A.size();
#ifdef DEBUG
    if(A.size() != B.size())
        throw std::runtime_error("mismatched sizes in MatRef/Mat call_daxpy");
    if(A.size() > std::numeric_limits<LAPACK_INT>::max()) 
        throw std::runtime_error("overflow of size beyond LAPACK_INT range");
#endif
    daxpy_wrapper(&size,&alpha,B.data(),&inc,A.data(),&inc);
    }

void
operator+=(const MatRef& a, MatRefc b)
    {
#ifdef DEBUG
    if(!(a.Ncols()==b.Ncols() && a.Nrows()==b.Nrows())) 
        throw std::runtime_error("MatRef +=: mismatched sizes");
#endif
    if(b.ind()==a.ind() && b.contiguous())
        {
        call_daxpy(a,b,+1);
        }
    else
        {
        auto pluseq = [](Real& x, Real y) { x += y; };
        apply(a,b.cbegin(),pluseq);
        }
    }

void
operator-=(const MatRef& a, MatRefc b)
    {
#ifdef DEBUG
    if(!(a.Ncols()==b.Ncols() && a.Nrows()==b.Nrows())) 
        throw std::runtime_error("MatRef +=: mismatched sizes");
#endif
    if(b.ind()==a.ind() && b.contiguous())
        {
        call_daxpy(a,b,-1);
        }
    else
        {
        auto minuseq = [](Real& x, Real y) { x -= y; };
        apply(a,b.cbegin(),minuseq);
        }
    }

Real
norm(MatRefc M)
    {
    Real nrm = 0;
    if(M.contiguous())
        {
        auto end = M.data()+M.size();
        for(auto it = M.data(); it != end; ++it) 
            {
            nrm += (*it)*(*it);
            }
        }
    else
        {
        for(auto& el : M) nrm += el*el;
        }
    return std::sqrt(nrm);
    }

void
randomize(MatRef M)
    {
    for(auto& el : M) el = detail::quickran();
    }


std::ostream&
operator<<(std::ostream& s, MatRefc M)
    {
    for(long r = 1; r <= M.Nrows(); ++r)
        {
        s << "|";
        for(long c = 1; c <= M.Ncols(); ++c)
            {
            s << M(r,c);
            s << (c == M.Ncols() ? "|" : " ");
            }
        if(r < M.Nrows()) s << "\n";
        }
    return s;
    }

// C = alpha*A*B + beta*C
void
call_dgemm(MatRefc A, 
           MatRefc B, 
           MatRef  C,
           Real alpha,
           Real beta)
    {
#ifdef DEBUG
    if(!(A.contiguous() && B.contiguous() && C.contiguous())) 
        throw std::runtime_error("multiplication of non-contiguous MatRefs not currently supported");
#endif
    if(C.transposed())
        {
        //Do C = Bt*At instead of Ct=A*B
        //Recall that C.data() points to elements of C, not C.t()
        //regardless of whether C.transpose()==true or false
        std::swap(A,B);
        A.applyTrans();
        B.applyTrans();
        C.applyTrans();
        }

#ifdef DEBUG
    if(A.Ncols() != B.Nrows())
        throw std::runtime_error("matrices A, B incompatible");
    if(A.Nrows() != C.Nrows() || B.Ncols() != C.Ncols())
        {
        printfln("A is %dx%d",A.Nrows(),A.Ncols());
        printfln("B is %dx%d",B.Nrows(),B.Ncols());
        printfln("C is %dx%d",C.Nrows(),C.Ncols());
        throw std::runtime_error("mult(_add) AxB -> C: matrix C incompatible");
        }
#endif
    dgemm_wrapper(A.transposed(),B.transposed(),
                  A.Nrows(),B.Ncols(),A.Ncols(),
                  alpha,A.data(),B.data(),beta,C.data());
    }

void
mult(MatRefc A, 
     MatRefc B, 
     MatRef  C)
    {
    call_dgemm(A,B,C,1.,0.);
    }

void
multAdd(MatRefc A, 
        MatRefc B, 
        MatRef  C)
    {
    call_dgemm(A,B,C,1.,1.);
    }

void
call_dgemv(const MatRefc& M,
           const VecRefc& x, 
           VecRef& y,
           Real alpha,
           Real beta,
           bool fromleft)
    {
#ifdef DEBUG
    if(!M.contiguous())
        throw std::runtime_error("multiplication of non-contiguous matrixref by vector not currently supported");
#endif
    auto trans = fromleft;
    LAPACK_INT m = M.Nrows(),
               n = M.Ncols();
    if(M.transposed())
        {
        trans = !fromleft;
        m = M.Ncols();
        n = M.Nrows();
        }
    dgemv_wrapper(trans,alpha,beta,m,n,M.data(),x.data(),x.stride(),y.data(),y.stride());
    }

void
mult(MatRefc M,
     VecRefc x,
     VecRef y,
     bool fromleft)
    {
#ifdef DEBUG
    if(fromleft ? M.Nrows()!=x.size() : M.Ncols()!=x.size()) 
        throw std::runtime_error("matrix vector mult: mismatched sizes");
    if(fromleft ? M.Ncols()!=y.size() : M.Nrows()!=y.size())
        throw std::runtime_error("matrix vector mult: wrong size for result (y) vec");
#endif
    call_dgemv(M,x,y,1,0,fromleft);
    }


} //namespace itensor
