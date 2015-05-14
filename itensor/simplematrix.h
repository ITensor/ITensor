//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SIMPLEMATRIX_H
#define __ITENSOR_SIMPLEMATRIX_H

#include "types.h"
#include "matrixref.h"
#include "print.h"

namespace itensor {

class SimpleMatrixRef
    {
    Real *store_ = nullptr;
    const Real *cstore_ = nullptr;
    long nrows_ = 0, 
         ncols_ = 0, 
         rowstride_ = 0;
    bool transpose_ = false;
    public:

    SimpleMatrixRef() { }

    SimpleMatrixRef(const Real* sto, 
                    long nro, 
                    long ncol, 
                    bool trans = false,
                    long rowstr = -1) 
        : 
        store_(nullptr), 
        cstore_(sto), 
        nrows_(nro), 
        ncols_(ncol), 
        rowstride_(rowstr < 0 ? ncol : rowstr), 
        transpose_(trans)
        { }

    SimpleMatrixRef(Real* sto, 
                    long nro, 
                    long ncol, 
                    bool trans = false,
                    long rowstr = -1) 
        : 
        SimpleMatrixRef(static_cast<const Real*>(sto),nro,ncol,trans,rowstr)
        { 
        store_ = sto;
        }

    SimpleMatrixRef(const SimpleMatrixRef& other) = default;

    SimpleMatrixRef(MatrixRefNoLink m)
        :
        store_(m.Store()), 
        cstore_(m.Store()), 
        nrows_(m.Nrows()), 
        ncols_(m.Ncols()), 
        rowstride_(m.RowStride()), 
        transpose_(m.DoTranspose())
        { }

    bool
    readOnly() const { return !bool(store_); }

    long
    Nrows() const { return transpose_ ? ncols_ : nrows_; }
    long
    Ncols() const { return transpose_ ? nrows_ : ncols_; }
    long
    rowStride() const { return rowstride_; }

    bool
    transpose() const { return transpose_; }
    void
    ApplyTrans() { transpose_ = !transpose_; }

    const Real*
    store() const { return cstore_; }
    void
    store(const Real* newstore) 
        { 
        store_ = nullptr;
        cstore_ = newstore; 
        }

    SimpleMatrixRef 
    t()
        { 
        SimpleMatrixRef res(*this);
        res.transpose_ = !transpose_; 
        return res;
        }

    Real
    operator()(long i, long j) const { return cstore_[index(i,j)]; }

    void
    set(Real val, long i, long j) 
        { 
#ifdef DEBUG
        if(readOnly()) Error("set: SimpleMatrixRef is read only");
#endif
        store_[index(i,j)] = val; 
        }

    private:
    long 
    index0(long i, long j) const
        { 
        return transpose_ ? j * rowstride_ + i : i * rowstride_ + j; 
        }
    long 
    index(long i, long j) const { return index0(i-1,j-1); }
    };

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
void 
mult_add(SimpleMatrixRef A, 
         SimpleMatrixRef B, 
         SimpleMatrixRef C, 
         Real beta = 1.0, 
         Real alpha = 1.0)
    {
#ifdef MATRIXBOUNDS
    if(A.Ncols() != B.Nrows())
        {
        printfln("A(%d,%d) * B(%d,%d) != C(%d,%d)",A.Nrows(),A.Ncols(),B.Nrows(),B.Ncols(),C.Nrows(),C.Ncols());
        Error("mult_add(A,B,C): Matrices A, B incompatible");
        }
    if(A.Nrows() != C.Nrows() || B.Ncols() != C.Ncols())
        {
        printfln("A(%d,%d) * B(%d,%d) != C(%d,%d)",A.Nrows(),A.Ncols(),B.Nrows(),B.Ncols(),C.Nrows(),C.Ncols());
        Error("mult_add(A,B,C): Matrix C incompatible");
        }
    if(C.readOnly()) Error("mult_add: error, C is readOnly");
#endif

    // Use BLAS 3 routine
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

    char transb = A.transpose() ? 'T' : 'N';
    char transa = B.transpose() ? 'T' : 'N';
    dgemm_(&transa,&transb,&m,&n,&k,&alpha,pa,&lda,pb,&ldb,&beta,pc,&ldc);
    }

} //namespace itensor

#endif
