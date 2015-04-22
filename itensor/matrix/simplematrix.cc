//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "simplematrix.h"

namespace itensor {

void vecref::
operator=(const vecref& other)
    {
    store_ = other.store_;
    cstore_ = other.cstore_;
    size_ = other.size_;
    }

matrixref::
matrixref(long nro, 
          long ncol, 
          bool trans)
    : 
    store_(nullptr),
    cstore_(nullptr),
    nrows_(nro), 
    ncols_(ncol), 
    transpose_(trans)
    { 
    }

matrixref::
matrixref(const Real* sto, 
          long nro, 
          long ncol, 
          bool trans)
    : 
    store_(nullptr),
    cstore_(sto),
    nrows_(nro), 
    ncols_(ncol), 
    transpose_(trans)
    { 
    }

matrixref::
matrixref(Real* sto, 
          long nro, 
          long ncol, 
          bool trans)
    : 
    store_(sto),
    cstore_(sto),
    nrows_(nro), 
    ncols_(ncol), 
    transpose_(trans)
    { 
    }

void matrixref::
operator=(const matrixref& other)
    {
    store_ = other.store_;
    cstore_ = other.cstore_;
    nrows_ = other.nrows_;
    ncols_ = other.ncols_;
    transpose_ = other.transpose_;
    }


matrixref matrixref::
t()
    { 
    matrixref res(*this);
    res.applyTrans();
    return res;
    }

std::ostream&
operator<<(std::ostream& s, const vecref& v)
    {
    for(auto j = 0ul; j <= v.size(); ++j)
        {
        s << v(j) << " ";
        }
    return s;
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

using BlasInt = int;
extern "C" void dgemm_(char*,char*,BlasInt*,BlasInt*,BlasInt*,Real*,Real*,BlasInt*,
	                   Real*,BlasInt*,Real*,Real*,BlasInt*);

// C = alpha*A*B + beta*C
void
dgemm_wrapper(const matrixref& AA, 
              const matrixref& BB, 
              matrixref& CC,
              Real beta,
              Real alpha)
    {
    //Workaround while getting slicing going:
    auto A = AA;
    auto B = BB;
    auto C = CC;

    if(C.transpose())
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
    BlasInt m = A.Nrows();
    BlasInt n = B.Ncols();
    BlasInt k = A.Ncols();
    char transa = A.transpose() ? 'T' : 'N';
    char transb = B.transpose() ? 'T' : 'N';
    BlasInt lda = A.transpose() ? k : m;
    BlasInt ldb = B.transpose() ? n : k;

    //New version:
    //BlasInt m = A.Nrows();
    //BlasInt n = B.Ncols();
    //BlasInt k = A.Ncols();
    //Something like:
    //char transa = A.index().rowMajor() ? 'N' : 'T';
    //char transb = ... similar
    //BlasInt lda = A.index().rowStride();
    //BlasInt ldb = B.index().rowStride();

    auto *pa = const_cast<Real*>(A.store());
    auto *pb = const_cast<Real*>(B.store());
    auto *pc = const_cast<Real*>(C.store());

    dgemm_(&transa,&transb,&m,&n,&k,&alpha,pa,&lda,pb,&ldb,&beta,pc,&m);
    }

void
mult_add(const matrixref& A, 
         const matrixref& B, 
         matrixref& C)
    {
    dgemm_wrapper(A,B,C,1,1);
    }

void
mult(const matrixref& A, 
     const matrixref& B, 
     matrixref& C)
    {
    dgemm_wrapper(A,B,C,0,1);
    }

}; //namespace itensor
