//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SIMPLEMATRIX_H_
#define __ITENSOR_SIMPLEMATRIX_H_

#include "types.h"
//#include "global.h"

namespace itensor {

class SimpleMatrixRef
    {
    const Real *store_ = nullptr;
    long nrows_ = 0, 
         ncols_ = 0;
    bool transpose_ = false;
    public:

    SimpleMatrixRef() { }

    SimpleMatrixRef(const Real* sto, 
                    long nro, 
                    long ncol, 
                    bool trans)
        : 
        store_(sto), 
        nrows_(nro), 
        ncols_(ncol), 
        transpose_(trans)
        { }

    SimpleMatrixRef(const SimpleMatrixRef& other) = default;

    long
    Nrows() const { return transpose_ ? ncols_ : nrows_; }
    long
    Ncols() const { return transpose_ ? nrows_ : ncols_; }
    long
    rowStride() const { return ncols_; }

    bool
    transpose() const { return transpose_; }
    void
    ApplyTrans() { transpose_ = !transpose_; }

    const Real*
    store() const { return store_; }
    void
    store(const Real* newstore) { store_ = newstore; }

    SimpleMatrixRef 
    t()
        { 
        SimpleMatrixRef res(*this);
        res.ApplyTrans();
        return res;
        }
    };

inline
std::ostream&
operator<<(std::ostream& s, const SimpleMatrixRef& M)
    {
    auto p = M.store();
    for(int r = 1; r <= M.Nrows(); ++r)
        {
        s << "|";
        for(int c = 1; c <= M.Ncols(); ++c)
            {
            s << (*p);
            s << (c == M.Ncols() ? "|" : " ");
            ++p;
            }
        s << "\n";
        }
    return s;
    }

using BlasInt = int;
extern "C" void dgemm_(char*,char*,BlasInt*,BlasInt*,BlasInt*,Real*,Real*,BlasInt*,
	                   Real*,BlasInt*,Real*,Real*,BlasInt*);

// C = alpha*A*B + beta*C
void inline
dgemm_wrapper(SimpleMatrixRef A, 
              SimpleMatrixRef B, 
              SimpleMatrixRef C, 
              Real beta,
              Real alpha)
    {
    // Use BLAS 3 routine
    if(!C.transpose())
        {
#ifdef DEBUG
        if(A.Ncols() != B.Nrows())
            throw std::runtime_error("mult_add(A,B,C): Matrices A, B incompatible");
        if(A.Nrows() != C.Nrows() || B.Ncols() != C.Ncols())
            throw std::runtime_error("mult_add(A,B,C): Matrix C incompatible");
#endif
        // Have to reverse the order, since we are really multiplying Ct = Bt*At
        BlasInt m = C.Ncols();
        BlasInt n = C.Nrows();
        BlasInt k = B.Nrows();
        BlasInt lda = B.rowStride();
        BlasInt ldb = A.rowStride();
        BlasInt ldc = C.rowStride();

        Real *pa = const_cast<Real*>(B.store());
        Real *pb = const_cast<Real*>(A.store());
        Real *pc = const_cast<Real*>(C.store());

        //printfln("Calling C no transpose version: m=%d, n=%d, k=%d, lda=%d, ldb=%d, B.t=%s, A.t=%s",m,n,k,lda,ldb,B.transpose(),A.transpose());
        char transa = B.transpose() ? 'T' : 'N';
        char transb = A.transpose() ? 'T' : 'N';
        dgemm_(&transa,&transb,&m,&n,&k,&alpha,pa,&lda,pb,&ldb,&beta,pc,&ldc);
        }
    else //C.transpose()==true
        {
#ifdef DEBUG
        if(A.Ncols() != B.Nrows())
            throw std::runtime_error("mult_add(A,B,C): Matrices A, B incompatible");
        if(A.Nrows() != C.Nrows() || B.Ncols() != C.Ncols())
            throw std::runtime_error("mult_add(A,B,C): Matrix C incompatible");
#endif
        // Here we are really multiplying C = A*B
        BlasInt m = C.Nrows();
        BlasInt n = C.Ncols();
        BlasInt k = A.Ncols();
        BlasInt lda = A.rowStride();
        BlasInt ldb = B.rowStride();
        BlasInt ldc = C.rowStride();

        Real *pa = const_cast<Real*>(A.store());
        Real *pb = const_cast<Real*>(B.store());
        Real *pc = const_cast<Real*>(C.store());

        //printfln("Calling C transpose version: m=%d, n=%d, k=%d, lda=%d, ldb=%d",m,n,k,lda,ldb);
        char transa = A.transpose() ? 'N' : 'T';
        char transb = B.transpose() ? 'N' : 'T';
        dgemm_(&transa,&transb,&m,&n,&k,&alpha,pa,&lda,pb,&ldb,&beta,pc,&ldc);
        }
    }

void inline
mult_add(SimpleMatrixRef A, 
         SimpleMatrixRef B, 
         SimpleMatrixRef C)
    {
    dgemm_wrapper(A,B,C,1,1);
    }

void inline
mult(SimpleMatrixRef A, 
     SimpleMatrixRef B, 
     SimpleMatrixRef C)
    {
    dgemm_wrapper(A,B,C,0,1);
    }

};

#endif
