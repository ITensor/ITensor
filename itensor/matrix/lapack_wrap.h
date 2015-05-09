//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_lapack_wrap_h
#define __ITENSOR_lapack_wrap_h

#include <vector>

//
// Headers and typedefs
//

#ifdef PLATFORM_lapack

#define LAPACK_REQUIRE_EXTERN

#include <complex>

namespace itensor {
using LAPACK_INT = int;
using LAPACK_REAL = double;
typedef struct
{
  double real, imag;
} LAPACK_COMPLEX;
};

#elif defined PLATFORM_macos

#include <Accelerate/Accelerate.h>
namespace itensor {
using LAPACK_INT = __CLPK_integer;
using LAPACK_REAL = __CLPK_doublereal;
using LAPACK_COMPLEX = __CLPK_doublecomplex;
};

#elif defined PLATFORM_acml

#define LAPACK_REQUIRE_EXTERN

//#include "acml.h"
namespace itensor {
using LAPACK_INT = int;
using LAPACK_REAL = double;
typedef struct
{
  double real, imag;
} LAPACK_COMPLEX;
};

#elif defined PLATFORM_mkl

#define FORTRAN_NO_TRAILING_UNDERSCORE

#include "mkl_blas.h"
#include "mkl_lapack.h"
namespace itensor {
using LAPACK_INT = MKL_INT;
using LAPACK_REAL = double;
using LAPACK_COMPLEX = MKL_Complex16;
};

#endif

#ifdef FORTRAN_NO_TRAILING_UNDERSCORE
#define F77NAME(x) x
#else
#define F77NAME(x) x##_
#endif

namespace itensor {

//
//
// Forward declarations of fortran lapack routines
//
//
#ifdef LAPACK_REQUIRE_EXTERN
extern "C" {


#ifdef PLATFORM_macos
void cblas_daxpy(const int n, const double alpha, const double *X, const int incX, double *Y, const int incY);
#else
void F77NAME(daxpy)(LAPACK_INT* n, LAPACK_REAL* alpha, 
                    LAPACK_REAL* X, LAPACK_INT* incx,
                    LAPACK_REAL* Y, LAPACK_INT* incy);
#endif

#ifdef PLATFORM_macos
LAPACK_REAL 
cblas_ddot(const LAPACK_INT N, const LAPACK_REAL *X, const LAPACK_INT incx, const LAPACK_REAL *Y, const LAPACK_INT incy);
#else
LAPACK_REAL F77NAME(ddot)(LAPACK_INT* N, LAPACK_REAL* X, LAPACK_INT* incx, LAPACK_REAL* Y, LAPACK_INT* incy);
#endif

#ifdef PLATFORM_macos
void cblas_dgemm(const enum CBLAS_ORDER __Order,
        const enum CBLAS_TRANSPOSE __TransA,
        const enum CBLAS_TRANSPOSE __TransB, const int __M, const int __N,
        const int __K, const double __alpha, const double *__A,
        const int __lda, const double *__B, const int __ldb,
        const double __beta, double *__C, const int __ldc);
#else
void F77NAME(dgemm)(char*,char*,LAPACK_INT*,LAPACK_INT*,LAPACK_INT*,
            LAPACK_REAL*,LAPACK_REAL*,LAPACK_INT*,LAPACK_REAL*,
            LAPACK_INT*,LAPACK_REAL*,LAPACK_REAL*,LAPACK_INT*);
#endif

#ifdef PLATFORM_macos
void cblas_dgemv(const enum CBLAS_ORDER Order,
        const enum CBLAS_TRANSPOSE TransA, const LAPACK_INT M, const LAPACK_INT N,
        const LAPACK_REAL alpha, const LAPACK_REAL *A, const LAPACK_INT lda,
        const LAPACK_REAL *X, const LAPACK_INT incX, const LAPACK_REAL beta, LAPACK_REAL *Y,
        const LAPACK_INT incY);
#else
void F77NAME(dgemv)(char* transa,LAPACK_INT* M,LAPACK_INT* N,LAPACK_REAL* alpha, LAPACK_REAL* A,
                    LAPACK_INT* LDA, LAPACK_REAL* X, LAPACK_INT* incx, LAPACK_REAL* beta,
                    LAPACK_REAL* Y, LAPACK_INT* incy);
#endif


#ifdef PLATFORM_acml
void F77NAME(dsyev)(char *jobz, char *uplo, int *n, double *a, int *lda, 
                    double *w, double *work, int *lwork, int *info, 
                    int jobz_len, int uplo_len);
#else
void F77NAME(dsyev)(const char* jobz, const char* uplo, const LAPACK_INT* n, double* a,
            const LAPACK_INT* lda, double* w, double* work, const LAPACK_INT* lwork,
            LAPACK_INT* info );
#endif
PA
#ifdef PLATFORM_macos
void cblas_dscal(const LAPACK_INT N, const LAPACK_REAL alpha, LAPACK_REAL* X,const LAPACK_INT incX);
#else
void F77NAME(dscal)(LAPACK_INT* N, LAPACK_REAL* alpha, LAPACK_REAL* X,LAPACK_INT* incX);
#endif

#ifdef PLATFORM_acml
void F77NAME(zgesdd)(char *jobz, int *m, int *n, LAPACK_COMPLEX *a, int *lda, double *s, 
             LAPACK_COMPLEX *u, int *ldu, LAPACK_COMPLEX *vt, int *ldvt, 
             LAPACK_COMPLEX *work, int *lwork, double *rwork, int *iwork, int *info, 
             int jobz_len);
#else
void F77NAME(zgesdd)(char *jobz, LAPACK_INT *m, LAPACK_INT *n, LAPACK_COMPLEX *a, LAPACK_INT *lda, double *s, 
             LAPACK_COMPLEX *u, LAPACK_INT *ldu, LAPACK_COMPLEX *vt, LAPACK_INT *ldvt, 
             LAPACK_COMPLEX *work, LAPACK_INT *lwork, double *rwork, LAPACK_INT *iwork, LAPACK_INT *info);
#endif

void F77NAME(dgeqrf)(LAPACK_INT *m, LAPACK_INT *n, double *a, LAPACK_INT *lda, 
                     double *tau, double *work, LAPACK_INT *lwork, LAPACK_INT *info);

void F77NAME(dorgqr)(LAPACK_INT *m, LAPACK_INT *n, LAPACK_INT *k, double *a, 
                     LAPACK_INT *lda, double *tau, double *work, LAPACK_INT *lwork, 
                     LAPACK_INT *info);

#ifdef PLATFORM_acml
void F77NAME(zheev)(char *jobz, char *uplo, LAPACK_INT *n, LAPACK_COMPLEX *a, LAPACK_INT *lda, 
            double *w, LAPACK_COMPLEX *work, LAPACK_INT *lwork, double *rwork, 
            LAPACK_INT *info, LAPACK_INT jobz_len, LAPACK_INT uplo_len);
#else
void F77NAME(zheev)(char *jobz, char *uplo, LAPACK_INT *n, LAPACK_COMPLEX *a, LAPACK_INT *lda, 
           double *w, LAPACK_COMPLEX *work, LAPACK_INT *lwork, double *rwork, 
           LAPACK_INT *info);
#endif


#ifdef PLATFORM_acml
void F77NAME(dsygv)(LAPACK_INT *itype, char *jobz, char *uplo, LAPACK_INT *n, double *a, 
            LAPACK_INT *lda, double *b, LAPACK_INT *ldb, double *w, double *work, 
            LAPACK_INT *lwork, LAPACK_INT *info, LAPACK_INT jobz_len, LAPACK_INT uplo_len); 
#else
void F77NAME(dsygv)(LAPACK_INT *itype, char *jobz, char *uplo, LAPACK_INT *n, double *a, 
           LAPACK_INT *lda, double *b, LAPACK_INT *ldb, double *w, double *work, 
           LAPACK_INT *lwork, LAPACK_INT *info);
#endif


#ifdef PLATFORM_acml
void F77NAME(dgeev)(char *jobvl, char *jobvr, LAPACK_INT *n, double *a, 
                    LAPACK_INT *lda, double *wr, double *wi, double *vl, LAPACK_INT *ldvl, 
                    double *vr, LAPACK_INT *ldvr, double *work, LAPACK_INT *lwork, 
                    LAPACK_INT *info, LAPACK_INT jobvl_len, LAPACK_INT jobvr_len);
#else
void F77NAME(dgeev)(char *jobvl, char *jobvr, LAPACK_INT *n, double *a, 
                    LAPACK_INT *lda, double *wr, double *wi, double *vl, LAPACK_INT *ldvl, 
                    double *vr, LAPACK_INT *ldvr, double *work, LAPACK_INT *lwork, 
                    LAPACK_INT *info);
#endif


#ifdef PLATFORM_acml
void F77NAME(zgeev)(char *jobvl, char *jobvr, LAPACK_INT *n, LAPACK_COMPLEX *a, 
                    LAPACK_INT *lda, LAPACK_COMPLEX *w, LAPACK_COMPLEX *vl, 
                    LAPACK_INT *ldvl, LAPACK_COMPLEX *vr, LAPACK_INT *ldvr, 
                    LAPACK_COMPLEX *work, LAPACK_INT *lwork, double *rwork, 
                    LAPACK_INT *info, LAPACK_INT jobvl_len, LAPACK_INT jobvr_len);
#else
void F77NAME(zgeev)(char *jobvl, char *jobvr, LAPACK_INT *n, LAPACK_COMPLEX *a, 
                    LAPACK_INT *lda, LAPACK_COMPLEX *w, LAPACK_COMPLEX *vl, 
                    LAPACK_INT *ldvl, LAPACK_COMPLEX *vr, LAPACK_INT *ldvr, 
                    LAPACK_COMPLEX *work, LAPACK_INT *lwork, double *rwork, 
                    LAPACK_INT *info);
#endif

} //extern "C"
#endif

//
// daxpy
// Y += alpha*X
//
void inline
daxpy_wrapper(LAPACK_INT* n,        //number of elements of X,Y
              LAPACK_REAL* alpha,   //scale factor
              const LAPACK_REAL* X, //pointer to head of vector X
              LAPACK_INT* incx,     //increment with which to step through X
              LAPACK_REAL* Y,       //pointer to head of vector Y
              LAPACK_INT* incy)     //increment with which to step through Y
    {
#ifdef PLATFORM_macos
    cblas_daxpy(*n,*alpha,X,*incx,Y,*incy);
#else
    auto Xnc = const_cast<LAPACK_REAL*>(X);
    F77NAME(daxpy)(n,alpha,Xnc,incx,Y,incy);
#endif
    }

//
// ddot
//
LAPACK_REAL inline
ddot_wrapper(LAPACK_INT N,
             const LAPACK_REAL* X,
             LAPACK_INT incx,
             const LAPACK_REAL* Y,
             LAPACK_INT incy)
    {
#ifdef PLATFORM_macos
    return cblas_ddot(N,X,incx,Y,incy);
#else
    auto *Xnc = const_cast<LAPACK_REAL*>(X);
    auto *Ync = const_cast<LAPACK_REAL*>(Y);
    return F77NAME(ddot)(&N,Xnc,&incx,Ync,&incy);
#endif
    return -1;
    }

//
// dgemm
//
void inline
dgemm_wrapper(bool transa, 
              bool transb,
              LAPACK_INT m,
              LAPACK_INT n,
              LAPACK_INT k,
              LAPACK_REAL alpha,
              const LAPACK_REAL* A,
              const LAPACK_REAL* B,
              LAPACK_REAL beta,
              LAPACK_REAL* C)
    {
    LAPACK_INT lda = m,
               ldb = k;
#ifdef PLATFORM_macos
    auto at = CblasNoTrans,
         bt = CblasNoTrans;
    if(transa)
        {
        at = CblasTrans;
        lda = k;
        }
    if(transb)
        {
        bt = CblasTrans;
        ldb = n;
        }
    cblas_dgemm(CblasColMajor,at,bt,m,n,k,alpha,A,lda,B,ldb,beta,C,m);
#else
    auto *pA = const_cast<double*>(A);
    auto *pB = const_cast<double*>(B);
    char at = 'N';
    char bt = 'N';
    if(transa)
        {
        at = 'T';
        lda = k;
        }
    if(transb)
        {
        bt = 'T';
        ldb = n;
        }
    F77NAME(dgemm)(&at,&bt,&m,&n,&k,&alpha,pA,&lda,pB,&ldb,&beta,C,&m);
#endif
    }

//
// dgemv - matrix*vector multiply
//
void inline
dgemv_wrapper(bool trans, 
              LAPACK_REAL alpha,
              LAPACK_REAL beta,
              LAPACK_INT m,
              LAPACK_INT n,
              const LAPACK_REAL* A,
              const LAPACK_REAL* x,
              LAPACK_INT incx,
              LAPACK_REAL* y,
              LAPACK_INT incy)
    {
#ifdef PLATFORM_macos
    auto Tr = trans ? CblasTrans : CblasNoTrans;
    cblas_dgemv(CblasColMajor,Tr,m,n,alpha,A,m,x,incx,beta,y,incy);
#else
    char Tr = trans ? 'T' : 'N';
    F77NAME(dgemv)(&Tr,&m,&n,&alpha,const_cast<LAPACK_REAL*>(A),&m,const_cast<LAPACK_REAL*>(x),&incx,&beta,y,&incy);
#endif
    }


//
// dsyev
//
void inline
dsyev_wrapper(char jobz,        //if jobz=='V', compute eigs and evecs
              char uplo,        //if uplo=='U', read from upper triangle of A
              LAPACK_INT n,     //number of cols of A
              LAPACK_REAL* A,    //symmetric matrix A
              LAPACK_REAL* eigs, //eigenvalues on return
              LAPACK_INT& info)  //error info
    {
    static std::vector<LAPACK_REAL> work;
    LAPACK_INT lda = n;

#ifdef PLATFORM_acml
    LAPACK_INT lwork = std::max(1,3*n-1);
    work.resize(lwork+2);
    F77NAME(dsyev)(&jobz,&uplo,&n,A,&lda,eigs,work.data(),&lwork,&info,1,1);
#else
    //Compute optimal workspace size (will be written to wkopt)
    LAPACK_INT lwork = -1; //tell dsyev to compute optimal size
    LAPACK_REAL wkopt = 0;
    F77NAME(dsyev)(&jobz,&uplo,&n,A,&lda,eigs,&wkopt,&lwork,&info);
    lwork = LAPACK_INT(wkopt);
    work.resize(lwork+2);
    F77NAME(dsyev)(&jobz,&uplo,&n,A,&lda,eigs,work.data(),&lwork,&info);
#endif
    }

//
// dscal
//
void inline
dscal_wrapper(LAPACK_INT N,
              LAPACK_REAL alpha,
              LAPACK_REAL* data,
              LAPACK_INT inc = 1)
    {
#ifdef PLATFORM_macos
    cblas_dscal(N,alpha,data,inc);
#else
    F77NAME(dscal)(&N,&alpha,data,&inc);
#endif
    }

void inline
zgesdd_wrapper(char *jobz,           //char* specifying how much of U, V to compute
                                     //choosing *jobz=='S' computes min(m,n) cols of U, V
               LAPACK_INT *m,        //number of rows of input matrix *A
               LAPACK_INT *n,        //number of cols of input matrix *A
               LAPACK_COMPLEX *A,    //contents of input matrix A
               LAPACK_REAL *s,       //on return, singular values of A
               LAPACK_COMPLEX *u,    //on return, unitary matrix U
               LAPACK_COMPLEX *vt,   //on return, unitary matrix V transpose
               LAPACK_INT *info)
    {
    static std::vector<LAPACK_COMPLEX> work;
    static std::vector<LAPACK_REAL> rwork;
    static std::vector<LAPACK_INT> iwork;
    LAPACK_INT l = std::min(*m,*n),
               g = std::max(*m,*n);
    LAPACK_INT lwork = l*l+2*l+g+100;
    work.resize(lwork);
    rwork.resize(5*l*(1+l));
    iwork.resize(8*l);
#ifdef PLATFORM_acml
    LAPACK_INT jobz_len = 1;
    F77NAME(zgesdd)(jobz,m,n,A,m,s,u,m,vt,n,work.data(),&lwork,rwork.data(),iwork.data(),info,jobz_len);
#else
    F77NAME(zgesdd)(jobz,m,n,A,m,s,u,m,vt,n,work.data(),&lwork,rwork.data(),iwork.data(),info);
#endif
    }

//
// dgeqrf
//
// QR factorization of a real matrix A
//
void inline
dgeqrf_wrapper(LAPACK_INT* m,     //number of rows of A
               LAPACK_INT* n,     //number of cols of A
               LAPACK_REAL* A,    //matrix A
                                  //on return upper triangle contains R
               LAPACK_INT* lda,   //size of A (usually same as n)
               LAPACK_REAL* tau,  //scalar factors of elementary reflectors
                                  //length should be min(m,n)
               LAPACK_INT* info)  //error info
    {
    static std::vector<LAPACK_REAL> work;
    int lwork = std::max(1,4*std::max(*n,*m));
    work.resize(lwork+2); 
    F77NAME(dgeqrf)(m,n,A,lda,tau,work.data(),&lwork,info);
    }

//
// dorgqr
//
// Generates Q from output of QR factorization routine dgeqrf (see above)
//
void inline
dorgqr_wrapper(LAPACK_INT* m,     //number of rows of A
               LAPACK_INT* n,     //number of cols of A
               LAPACK_INT* k,     //number of elementary reflectors, typically min(m,n)
               LAPACK_REAL* A,    //matrix A, as returned from "A" argument of dgeqrf
                                  //on return contains Q
               LAPACK_INT* lda,   //size of A (usually same as n)
               LAPACK_REAL* tau,  //scalar factors as returned by dgeqrf
               LAPACK_INT* info)  //error info
    {
    static std::vector<LAPACK_REAL> work;
    auto lwork = std::max(1,4*std::max(*n,*m));
    work.resize(lwork+2); 
    F77NAME(dorgqr)(m,n,k,A,lda,tau,work.data(),&lwork,info);
    }

//
// zheev
//
// Eigenvalues and eigenvectors of complex Hermitian matrix A
//
void inline
zheev_wrapper(char* jobz,           //if 'V', compute both eigs and evecs
                                    //if 'N', only eigenvalues
              char* uplo,           //if 'U', use upper triangle of A
              LAPACK_INT* n,        //number of cols of A
              LAPACK_COMPLEX* A,    //matrix A, on return contains eigenvectors
              LAPACK_INT* lda,      //size of A (usually same as n)
              LAPACK_REAL* d,       //eigenvalues on return
              LAPACK_COMPLEX* work, //complex workspace array
              LAPACK_INT* lwork,    //size of work
              LAPACK_REAL* rwork,   //real workspace array
              LAPACK_INT* info)  //error info
    {
#ifdef PLATFORM_acml
    LAPACK_INT jobz_len = 1;
    LAPACK_INT uplo_len = 1;
    F77NAME(zheev)(jobz,uplo,n,A,lda,d,work,lwork,rwork,info,jobz_len,uplo_len);
#else
    F77NAME(zheev)(jobz,uplo,n,A,lda,d,work,lwork,rwork,info);
#endif
    }

//
// dsygv
//
// Eigenvalues and eigenvectors of generalized eigenvalue problem
// A*x = lambda*B*x
// A and B must be symmetric
// B must be positive definite
//
void inline
dsygv_wrapper(char* jobz,           //if 'V', compute both eigs and evecs
                                    //if 'N', only eigenvalues
              char* uplo,           //if 'U', use upper triangle of A
              LAPACK_INT* n,        //number of cols of A
              LAPACK_REAL* A,       //matrix A, on return contains eigenvectors
              LAPACK_REAL* B,       //matrix B
              LAPACK_REAL* d,       //eigenvalues on return
              LAPACK_INT* info)  //error info
    {
    static std::vector<LAPACK_REAL> work;
    int itype = 1;
    LAPACK_INT lwork = std::max(1,3*(*n)-1);//std::max(1, 1+6*N+2*N*N);
    work.resize(lwork);
#ifdef PLATFORM_acml
    LAPACK_INT jobz_len = 1;
    LAPACK_INT uplo_len = 1;
    F77NAME(dsygv)(&itype,jobz,uplo,n,A,n,B,n,d,work.data(),&lwork,info,jobz_len,uplo_len);
#else
    F77NAME(dsygv)(&itype,jobz,uplo,n,A,n,B,n,d,work.data(),&lwork,info);
#endif
    }

//
// dgeev
//
// Eigenvalues and eigenvectors of real, square matrix A
// A can be a general real matrix, not assumed symmetric
//
void inline
dgeev_wrapper(char* jobvl,          //if 'V', compute left eigenvectors, else 'N'
              char* jobvr,          //if 'V', compute right eigenvectors, else 'N'
              LAPACK_INT* n,        //number of rows/cols of A
              LAPACK_REAL* A,       //matrix A, on return contains eigenvectors
              LAPACK_REAL* dr,      //real parts of eigenvalues
              LAPACK_REAL* di,      //imaginary parts of eigenvalues
              LAPACK_REAL* vl,      //left eigenvectors on return
              LAPACK_REAL* vr,      //right eigenvectors on return
              LAPACK_INT* info)  //error info
    {
    static std::vector<LAPACK_REAL> work;
    LAPACK_INT nevecl = (*jobvl == 'V' ? *n : 1);
    LAPACK_INT nevecr = (*jobvr == 'V' ? *n : 1);
    LAPACK_INT lwork = std::max(1,4*(*n));
    work.resize(lwork);
#ifdef PLATFORM_acml
    LAPACK_INT jobvl_len = 1;
    LAPACK_INT jobvr_len = 1;
    F77NAME(dgeev)(jobvl,jobvr,n,A,n,dr,di,vl,&nevecl,vr,&nevecr,work.data(),&lwork,info,jobvl_len,jobvr_len);
#else
    F77NAME(dgeev)(jobvl,jobvr,n,A,n,dr,di,vl,&nevecl,vr,&nevecr,work.data(),&lwork,info);
#endif
    }

//
// zgeev
//
// Eigenvalues and eigenvectors of complex, square matrix A
// A can be a general complex matrix, not assumed symmetric
//
void inline
zgeev_wrapper(char* jobvl,          //if 'V', compute left eigenvectors, else 'N'
              char* jobvr,          //if 'V', compute right eigenvectors, else 'N'
              LAPACK_INT* n,        //number of rows/cols of A
              LAPACK_COMPLEX* A,    //matrix A, on return contains eigenvectors
              LAPACK_COMPLEX* d,    //eigenvalues
              LAPACK_COMPLEX* vl,   //left eigenvectors on return
              LAPACK_COMPLEX* vr,   //right eigenvectors on return
              LAPACK_INT* info)  //error info
    {
    static std::vector<LAPACK_COMPLEX> work;
    static std::vector<LAPACK_REAL> rwork;
    int nevecl = (*jobvl == 'V' ? *n : 1);
    int nevecr = (*jobvr == 'V' ? *n : 1);
    LAPACK_INT lwork = std::max(1,4*(*n));
    work.resize(lwork);
    LAPACK_INT lrwork = std::max(1,2*(*n));
    rwork.resize(lrwork);
#ifdef PLATFORM_acml
    F77NAME(zgeev)(jobvl,jobvr,n,A,n,d,vl,&nevecl,vr,&nevecr,work.data(),&lwork,rwork.data(),info,1,1);
#else
    F77NAME(zgeev)(jobvl,jobvr,n,A,n,d,vl,&nevecl,vr,&nevecr,work.data(),&lwork,rwork.data(),info);
#endif
    }

}; //namespace itensor

#endif
