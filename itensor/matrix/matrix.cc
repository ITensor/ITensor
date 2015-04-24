//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "matrix.h"
#include "lapack_wrap.h"
#include <limits>

namespace itensor {

//extern "C" void dgemm_(char*,char*,LAPACK_INT*,LAPACK_INT*,LAPACK_INT*,
//            LAPACK_REAL*,LAPACK_REAL*,LAPACK_INT*,LAPACK_REAL*,
//            LAPACK_INT*,LAPACK_REAL*,LAPACK_REAL*,LAPACK_INT*);

matrixref::
matrixref(long nro, 
          long ncol, 
          bool trans)
    : 
    ind_(trans ? mrange(ncol,nro,nro,1) : mrange(nro,ncol)),
    store_(nullptr),
    cstore_(nullptr)
    { 
    }

matrixref::
matrixref(const Real* sto, 
          long nro, 
          long ncol, 
          bool trans)
    : 
    ind_(trans ? mrange(ncol,nro,nro,1) : mrange(nro,ncol)),
    store_(nullptr),
    cstore_(sto)
    { 
    }

matrixref::
matrixref(Real* sto, 
          long nro, 
          long ncol, 
          bool trans)
    : 
    ind_(trans ? mrange(ncol,nro,nro,1) : mrange(nro,ncol)),
    store_(sto),
    cstore_(sto)
    { 
    }

matrixref::
matrixref(const Real* sto, 
          const mrange& ind)
    : 
    ind_(ind),
    store_(nullptr),
    cstore_(sto)
    { 
    }

matrixref::
matrixref(Real* sto, 
          const mrange& ind)
    : 
    ind_(ind),
    store_(sto),
    cstore_(sto)
    { 
    }

void matrixref::
operator=(const matrixref& other)
    {
    ind_ = other.ind_;
    store_ = other.store_;
    cstore_ = other.cstore_;
    }


matrixref matrixref::
t()
    { 
    matrixref res(*this);
    res.applyTrans();
    return res;
    }

std::ostream&
operator<<(std::ostream& s, const matrixref& M)
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

vecref
diagonal(const matrixref& m) 
    { 
    auto vsize = std::min(m.Nrows(),m.Ncols());
    auto vstrd = m.rowStride()+m.colStride();
    if(m.readOnly()) return vecref(m.cstore(),vsize,vstrd);
    else             return vecref(m.store(),vsize,vstrd);
    }

matrixref
subMatrix(const matrixref& m,
          long rstart,
          long rstop,
          long cstart,
          long cstop)
    { 
#ifdef DEBUG
    if(rstop > m.Nrows() || rstart >= rstop) throw std::runtime_error("subMatrix invalid row start and stop");
    if(cstop > m.Ncols() || cstart >= cstop) throw std::runtime_error("subMatrix invalid col start and stop");
#endif
    const auto& i = m.ind();
    auto offset = i.rs*(rstart-1)+i.cs*(cstart-1);
    auto subind = mrange(rstop-rstart+1,i.rs,cstop-cstart+1,i.cs);
    if(m.readOnly()) return matrixref(m.cstore()+offset,subind);
    else             return matrixref(m.store()+offset,subind);
    }


// C = alpha*A*B + beta*C
void
call_dgemm(matrixref A, 
           matrixref B, 
           matrixref C,
           Real beta,
           Real alpha)
    {
#ifdef DEBUG
    if(C.readOnly()) throw std::runtime_error("Can't store result of matrix multiply in read-only matrixref or matrix");
    if(!(A.contiguous() && B.contiguous() && C.contiguous())) 
        throw std::runtime_error("multiplication of non-contiguous matrixref's not currently supported");
#endif

    if(C.transposed())
        {
        //Do C = Bt*At instead of Ct=A*B
        //Recall that C.store() points to elements of C, not C.t()
        //regardless of whether C.transpose()==true or false
        std::swap(A,B);
        A.applyTrans();
        B.applyTrans();
        C.applyTrans();
        }

#ifdef DEBUG
    if(A.Ncols() != B.Nrows())
        throw std::runtime_error("mult_add(A,B,C): Matrices A, B incompatible");
    if(A.Nrows() != C.Nrows() || B.Ncols() != C.Ncols())
        {
        printfln("A is %dx%d",A.Nrows(),A.Ncols());
        printfln("B is %dx%d",B.Nrows(),B.Ncols());
        printfln("C is %dx%d",C.Nrows(),C.Ncols());
        throw std::runtime_error("mult_add(A,B,C): Matrix C incompatible");
        }
#endif
    LAPACK_INT m = A.Nrows();
    LAPACK_INT n = B.Ncols();
    LAPACK_INT k = A.Ncols();
    auto at = A.transposed();
    LAPACK_INT lda = at ? A.rowStride() : A.colStride();
    auto bt = B.transposed();
    LAPACK_INT ldb = bt ? B.rowStride() : B.colStride();

    auto *pa = const_cast<Real*>(A.cstore());
    auto *pb = const_cast<Real*>(B.cstore());
    auto *pc = C.store();

    dgemm_wrapper(at,bt,m,n,k,alpha,pa,lda,pb,ldb,beta,pc,m);
    }

void
mult_add(const matrixref& A, 
         const matrixref& B, 
         matrixref& C)
    {
    call_dgemm(A,B,C,1,1);
    }

void
mult(const matrixref& A, 
     const matrixref& B, 
     matrixref& C)
    {
    call_dgemm(A,B,C,0,1);
    }


matrix
operator+(const matrixref& A, const matrixref& B)
    {
#ifdef DEBUG
    if(!(A.Nrows() == B.Nrows() && A.Ncols() == B.Ncols())) throw std::runtime_error("Mismatched matrix sizes in plus");
#endif
    matrix C(A.Nrows(),A.Ncols());
//    if(A.ind() == B.ind())
//        {
//        LAPACK_REAL alpha = 1.0;
//        LAPACK_INT inc = 1;
//        LAPACK_INT size = C.size();
//#ifdef DEBUG
//        if(C.size() > std::numeric_limits<LAPACK_INT>::max()) throw std::runtime_error("plus: overflow of size beyond LAPACK_INT range");
//#endif
//        daxpy_wrapper(&size,&alpha,B.cstore(),&inc,C.store(),&inc);
//        return C;
//        }
    auto c = C.begin();
    auto a = A.cbegin();
    for(const auto& el : B)
        {
        *c = el + *a;
        ++c;
        ++a;
        }
    return C;
    }

void
diagSymmetric(const matrixref& M,
              matrixref& U,
              vecref& d)

    {
    LAPACK_INT N = M.Ncols();
    if(N < 1) throw std::runtime_error("diagSymmetric: 0 dimensional matrix");
    if(N != M.Nrows())
        {
        printfln("M is %dx%d",M.Nrows(),M.Ncols());
        throw std::runtime_error("diagSymmetric: Input Matrix must be square");
        }

#ifdef DEBUG
    if(U.Nrows() != N || U.Ncols() != N) 
        throw std::runtime_error("diagSymmetric: U should have same dims as M");
    if(!U.contiguous())
        throw std::runtime_error("diagSymmetric: U must be contiguous");
    if(!d.contiguous())
        throw std::runtime_error("diagSymmetric: d must be contiguous");
#endif

    char jobz = 'V';
    char uplo = 'U';
    LAPACK_INT info;

    std::copy(M.cbegin(),M.cend(),U.begin());
    
    dsyev_wrapper(&jobz,&uplo,&N,U.store(),&N,d.store(),&info);

    if(info != 0) throw std::runtime_error("Error condition in diagSymmetric");

    //Transpose U before return
    U.applyTrans();
    }

void
diagSymmetric(const matrixref& M,
              matrix& U,
              vec& d)
    {
    if(U.Nrows() != M.Nrows() || U.Ncols() != M.Ncols())
        U = matrix(M.Nrows(),M.Ncols());
    if(d.size() != M.Nrows()) d = vec(M.Nrows());
    matrixref& Uref = U;
    vecref& dref = d;
    diagSymmetric(M,Uref,dref);
    }

}; //namespace itensor
