//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "matrix.h"
#include "lapack_wrap.h"
#include <limits>

namespace itensor {

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
call_dgemm(const matrixref& AA, 
           const matrixref& BB, 
           matrixref& CC,
           Real beta,
           Real alpha)
    {
    //Need to modify A,B,C below but want to take
    //advantage of subtle C++ feature where const references
    //extend lifetimes of temporaries e.g. if AA is a temporary matrix object
    auto A = AA;
    auto B = BB;
    auto C = CC;
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


//matrix
//operator+(const matrixref& A, const matrixref& B)
//    {
//#ifdef DEBUG
//    if(!(A.Nrows() == B.Nrows() && A.Ncols() == B.Ncols())) throw std::runtime_error("Mismatched matrix sizes in plus");
//#endif
//    matrix C(A.Nrows(),A.Ncols());
//    auto c = C.begin();
//    auto a = A.cbegin();
//    for(const auto& el : B)
//        {
//        *c = el + *a;
//        ++c;
//        ++a;
//        }
//    return C;
//    }

void matrix::
call_daxpy(const matrix& other, Real alpha_)
    {
    LAPACK_REAL alpha = alpha_;
    LAPACK_INT inc = 1;
    LAPACK_INT size = parent::size();
#ifdef DEBUG
    if(!(other.Nrows() == Nrows() && other.Ncols() == Ncols())) throw std::runtime_error("matrix +=: mismatches sizes");
    if(parent::size() > std::numeric_limits<LAPACK_INT>::max()) throw std::runtime_error("matrix +=: overflow of size beyond LAPACK_INT range");
#endif
    daxpy_wrapper(&size,&alpha,other.cstore(),&inc,store(),&inc);
    }

matrix& matrix::
operator*=(Real fac)
    {
    LAPACK_INT sz = parent::size();
#ifdef DEBUG
    if(parent::size() > std::numeric_limits<LAPACK_INT>::max()) 
        throw std::runtime_error("matrix *=: overflow of size beyond LAPACK_INT range");
#endif
    dscal_wrapper(sz,fac,store());
    return *this;
    }
matrix& matrix::
operator/=(Real fac)
    {
    if(fac == 0) throw std::runtime_error("matrix /=: divide by zero");
    return operator*=(1./fac);
    }

}; //namespace itensor
