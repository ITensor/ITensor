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
#ifndef __ITENSOR_LAPACK_WRAP_h
#define __ITENSOR_LAPACK_WRAP_h

#include <vector>
#include "itensor/config.h"
#include "itensor/types.h"
#include "itensor/util/timers.h"

//
// Headers and typedefs
//

//
//
// Generic Linux LAPACK
//
//
#ifdef PLATFORM_lapack

#define LAPACK_REQUIRE_EXTERN

namespace itensor {
    using LAPACK_INT = int;
    using LAPACK_REAL = double;
    typedef struct
    {
    LAPACK_REAL real, imag;
    } LAPACK_COMPLEX;
}
#elif defined PLATFORM_openblas

#define ITENSOR_USE_CBLAS

#include "cblas.h"
#include "lapacke.h"
#undef I //lapacke.h includes complex.h which defined an `I` macro
         //that can cause problems, so best to undefine it

namespace itensor {
using LAPACK_INT = lapack_int;
using LAPACK_REAL = double;
using LAPACK_COMPLEX = lapack_complex_double;

inline LAPACK_REAL& 
realRef(LAPACK_COMPLEX & z) 
    { 
    auto* p = reinterpret_cast<double*>(&z);
    return p[0];
    }

inline LAPACK_REAL& 
imagRef(LAPACK_COMPLEX & z) 
    { 
    auto* p = reinterpret_cast<double*>(&z);
    return p[1];
    }
}

//
//
// Apple Accelerate/vecLib
//
//
#elif defined PLATFORM_macos

#define ITENSOR_USE_CBLAS
//#define ITENSOR_USE_ZGEMM

#include <Accelerate/Accelerate.h>
    namespace itensor {
    using LAPACK_INT = __CLPK_integer;
    using LAPACK_REAL = __CLPK_doublereal;
    using LAPACK_COMPLEX = __CLPK_doublecomplex;

    inline LAPACK_REAL& 
    realRef(LAPACK_COMPLEX & z) { return z.r; }

    inline LAPACK_REAL& 
    imagRef(LAPACK_COMPLEX & z) { return z.i; }
    }

//
//
// Intel MKL
//
//
#elif defined PLATFORM_mkl

#define ITENSOR_USE_CBLAS
#define ITENSOR_USE_ZGEMM

#include "mkl_cblas.h"
#include "mkl_lapack.h"
    namespace itensor {
    using LAPACK_INT = MKL_INT;
    using LAPACK_REAL = double;
    using LAPACK_COMPLEX = MKL_Complex16;

    inline LAPACK_REAL& 
    realRef(LAPACK_COMPLEX & z) { return z.real; }

    inline LAPACK_REAL& 
    imagRef(LAPACK_COMPLEX & z) { return z.imag; }
    }

//
//
// AMD ACML
//
//
#elif defined PLATFORM_acml

#define LAPACK_REQUIRE_EXTERN
//#include "acml.h"
    namespace itensor {
    using LAPACK_INT = int;
    using LAPACK_REAL = double;
    typedef struct
    {
    LAPACK_REAL real, imag;
    } LAPACK_COMPLEX;

    inline LAPACK_REAL& 
    realRef(LAPACK_COMPLEX & z) { return z.real; }

    inline LAPACK_REAL& 
    imagRef(LAPACK_COMPLEX & z) { return z.imag; }
    }

#endif // different PLATFORM types



#ifdef FORTRAN_NO_TRAILING_UNDERSCORE
#define F77NAME(x) x
#else
#if defined(LAPACK_GLOBAL) || defined(LAPACK_NAME)
#define F77NAME(x) LAPACK_##x
#else
#define F77NAME(x) x##_
#endif
#endif

namespace itensor {

//
//
// Forward declarations of fortran lapack routines
//
//
#ifdef LAPACK_REQUIRE_EXTERN
extern "C" {

//dnrm2 declaration
#ifdef ITENSOR_USE_CBLAS
LAPACK_REAL cblas_dnrm2(LAPACK_INT N, LAPACK_REAL *X, LAPACK_INT incX);
#else
LAPACK_REAL F77NAME(dnrm2)(LAPACK_INT* N, LAPACK_REAL* X, LAPACK_INT* incx);
#endif


//daxpy declaration
#ifdef ITENSOR_USE_CBLAS
void cblas_daxpy(const int n, const double alpha, const double *X, const int incX, double *Y, const int incY);
#else
void F77NAME(daxpy)(LAPACK_INT* n, LAPACK_REAL* alpha, 
                    LAPACK_REAL* X, LAPACK_INT* incx,
                    LAPACK_REAL* Y, LAPACK_INT* incy);
#endif

//ddot declaration
#ifdef ITENSOR_USE_CBLAS
LAPACK_REAL 
cblas_ddot(const LAPACK_INT N, const LAPACK_REAL *X, const LAPACK_INT incx, const LAPACK_REAL *Y, const LAPACK_INT incy);
#else
LAPACK_REAL F77NAME(ddot)(LAPACK_INT* N, LAPACK_REAL* X, LAPACK_INT* incx, LAPACK_REAL* Y, LAPACK_INT* incy);
#endif

//zdotc declaration
#ifdef ITENSOR_USE_CBLAS
LAPACK_REAL 
cblas_zdotc_sub(const LAPACK_INT N, const void *X, const LAPACK_INT incx, const void *Y, const LAPACK_INT incy, void *res);
#else
LAPACK_COMPLEX F77NAME(zdotc)(LAPACK_INT* N, LAPACK_COMPLEX* X, LAPACK_INT* incx, LAPACK_COMPLEX* Y, LAPACK_INT* incy);
#endif

//dgemm declaration
#ifdef ITENSOR_USE_CBLAS
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

//zgemm declaration
#ifdef PLATFORM_openblas
void cblas_zgemm(OPENBLAS_CONST enum CBLAS_ORDER Order, 
                 OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, 
                 OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, 
                 OPENBLAS_CONST blasint M, 
                 OPENBLAS_CONST blasint N, 
                 OPENBLAS_CONST blasint K,
                 OPENBLAS_CONST double *alpha, 
                 OPENBLAS_CONST double *A, 
                 OPENBLAS_CONST blasint lda, 
                 OPENBLAS_CONST double *B, 
                 OPENBLAS_CONST blasint ldb, 
                 OPENBLAS_CONST double *beta, 
                 double *C, 
                 OPENBLAS_CONST blasint ldc);
#else //platform not openblas

#ifdef ITENSOR_USE_CBLAS
void cblas_zgemm(const enum CBLAS_ORDER __Order,
        const enum CBLAS_TRANSPOSE __TransA,
        const enum CBLAS_TRANSPOSE __TransB, const int __M, const int __N,
        const int __K, const void *__alpha, const void *__A, const int __lda,
        const void *__B, const int __ldb, const void *__beta, void *__C,
        const int __ldc);
#else
void F77NAME(zgemm)(char* transa,char* transb,LAPACK_INT* m,LAPACK_INT* n,LAPACK_INT* k,
            LAPACK_COMPLEX* alpha,LAPACK_COMPLEX* A,LAPACK_INT* LDA,LAPACK_COMPLEX* B,
            LAPACK_INT* LDB,LAPACK_COMPLEX* beta,LAPACK_COMPLEX* C,LAPACK_INT* LDC);
#endif

#endif //zgemm declaration

//dgemv declaration
#ifdef ITENSOR_USE_CBLAS
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

//zgemv declaration
#ifdef PLATFORM_openblas
void cblas_zgemv(OPENBLAS_CONST enum CBLAS_ORDER order,  
                 OPENBLAS_CONST enum CBLAS_TRANSPOSE trans,  
                 OPENBLAS_CONST blasint m, 
                 OPENBLAS_CONST blasint n,
                 OPENBLAS_CONST double *alpha, 
                 OPENBLAS_CONST double  *a, 
                 OPENBLAS_CONST blasint lda,  
                 OPENBLAS_CONST double  *x, 
                 OPENBLAS_CONST blasint incx,  
                 OPENBLAS_CONST double *beta,  
                 double  *y, 
                 OPENBLAS_CONST blasint incy);
#else
#ifdef ITENSOR_USE_CBLAS
void cblas_zgemv(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE trans, const LAPACK_INT m, 
                 const LAPACK_INT n, const void *alpha, const void *a, const LAPACK_INT lda, 
                 const void *x, const LAPACK_INT incx, const void *beta, void *y, const LAPACK_INT incy);
#else
void F77NAME(zgemv)(char* transa,LAPACK_INT* M,LAPACK_INT* N,LAPACK_COMPLEX* alpha, LAPACK_COMPLEX* A,
                    LAPACK_INT* LDA, LAPACK_COMPLEX* X, LAPACK_INT* incx, LAPACK_COMPLEX* beta,
                    LAPACK_COMPLEX* Y, LAPACK_INT* incy);
#endif
#endif //zgemv declaration

#ifdef PLATFORM_acml
void F77NAME(dsyev)(char *jobz, char *uplo, int *n, double *a, int *lda, 
                    double *w, double *work, int *lwork, int *info, 
                    int jobz_len, int uplo_len);
#else
void F77NAME(dsyev)(const char* jobz, const char* uplo, const LAPACK_INT* n, double* a,
            const LAPACK_INT* lda, double* w, double* work, const LAPACK_INT* lwork,
            LAPACK_INT* info );
#endif

#ifdef ITENSOR_USE_CBLAS
void cblas_dscal(const LAPACK_INT N, const LAPACK_REAL alpha, LAPACK_REAL* X,const LAPACK_INT incX);
#else
void F77NAME(dscal)(LAPACK_INT* N, LAPACK_REAL* alpha, LAPACK_REAL* X,LAPACK_INT* incX);
#endif


#ifdef PLATFORM_acml
void F77NAME(dgesdd)(char *jobz, LAPACK_INT *m, LAPACK_INT *n, double *a, LAPACK_INT *lda, double *s, 
             double *u, LAPACK_INT *ldu, double *vt, LAPACK_INT *ldvt, 
             double *work, LAPACK_INT *lwork, LAPACK_INT *iwork, LAPACK_INT *info, int jobz_len);
#else
void F77NAME(dgesdd)(char *jobz, LAPACK_INT *m, LAPACK_INT *n, double *a, LAPACK_INT *lda, double *s, 
             double *u, LAPACK_INT *ldu, double *vt, LAPACK_INT *ldvt, 
             double *work, LAPACK_INT *lwork, LAPACK_INT *iwork, LAPACK_INT *info);
#endif


#ifdef PLATFORM_acml
  void F77NAME(dgesvd)(char *jobz, char* jobv, LAPACK_INT *m, LAPACK_INT *n, double *a, LAPACK_INT *lda, double *s, 
             double *u, LAPACK_INT *ldu, double *vt, LAPACK_INT *ldvt, 
             double *work, LAPACK_INT *lwork, LAPACK_INT *info, int jobz_len);
#else
  void F77NAME(dgesvd)(char *jobz, char* jobv, LAPACK_INT *m, LAPACK_INT *n, double *a, LAPACK_INT *lda, double *s, 
             double *u, LAPACK_INT *ldu, double *vt, LAPACK_INT *ldvt, 
             double *work, LAPACK_INT *lwork, LAPACK_INT *info);
#endif


  #ifdef PLATFORM_acml
  void F77NAME(zgesvd)(char *jobz, char* jobv, LAPACK_INT *m, LAPACK_INT *n, LAPACK_COMPLEX *a, LAPACK_INT *lda, LAPACK_REAL *s, 
             LAPACK_COMPLEX *u, LAPACK_INT *ldu,  LAPACK_COMPLEX *vt, LAPACK_INT *ldvt, 
             LAPACK_COMPLEX *work, LAPACK_INT *lwork, LAPACK_REAL * rwork, LAPACK_INT *info, int jobz_len);
#else
  void F77NAME(zgesvd)(char *jobz, char* jobv, LAPACK_INT *m, LAPACK_INT *n, LAPACK_COMPLEX *a, LAPACK_INT *lda, LAPACK_REAL *s, 
             LAPACK_COMPLEX *u, LAPACK_INT *ldu, LAPACK_COMPLEX *vt, LAPACK_INT *ldvt, 
		       LAPACK_COMPLEX *work, LAPACK_INT *lwork, LAPACK_REAL * rwork, LAPACK_INT *info);
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

  
void F77NAME(zgeqrf)(LAPACK_INT *m, LAPACK_INT *n, LAPACK_COMPLEX *a, LAPACK_INT *lda, 
                     LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INT *lwork, LAPACK_INT *info);

#ifdef PLATFORM_lapacke
void LAPACKE_zungqr(int matrix_layout, LAPACK_INT *m, LAPACK_INT *n, LAPACK_INT *k, LAPACK_COMPLEX *a, 
                     LAPACK_INT *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INT *lwork, 
                     LAPACK_INT *info);
#else
void F77NAME(zungqr)(LAPACK_INT *m, LAPACK_INT *n, LAPACK_INT *k, LAPACK_COMPLEX *a, 
                     LAPACK_INT *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INT *lwork, 
                     LAPACK_INT *info);
#endif

void F77NAME(dgesv)(LAPACK_INT *n, LAPACK_INT *nrhs, LAPACK_REAL *a, LAPACK_INT *lda,
					LAPACK_INT *ipiv, LAPACK_REAL *b, LAPACK_INT *ldb, LAPACK_INT *info);

void F77NAME(zgesv)(LAPACK_INT *n, LAPACK_INT *nrhs, LAPACK_COMPLEX *a, LAPACK_INT * lda,
					LAPACK_INT *ipiv, LAPACK_COMPLEX *b, LAPACK_INT *ldb, LAPACK_INT *info);

#ifdef PLATFORM_lapacke
double LAPACKE_dlange(int matrix_layout, char norm, lapack_int m, lapack_int n, const double* a, lapack_int lda);
#elif defined PLATFORM_acml
double F77NAME(dlange)(char* norm, LAPACK_INT* m, LAPACK_INT* n, double* a, LAPACK_INT* lda, double* work, LAPACK_INT norm_len);
#else
double F77NAME(dlange)(char* norm, LAPACK_INT* m, LAPACK_INT* n, double* a, LAPACK_INT* lda, double* work);
#endif

#ifdef PLATFORM_lapacke
lapack_real LAPACKE_zlange(int matrix_layout, char norm, lapack_int m, lapack_int n, const lapack_complex_double* a, lapack_int lda);
#elif defined PLATFORM_acml
LAPACK_REAL F77NAME(zlange)(char* norm, LAPACK_INT* m, LAPACK_INT* n, LAPACK_COMPLEX* a, LAPACK_INT* lda, double* work, LAPACK_INT norm_len);
#else
LAPACK_REAL F77NAME(zlange)(char* norm, LAPACK_INT* m, LAPACK_INT* n, LAPACK_COMPLEX* a, LAPACK_INT* lda, double* work);
#endif

#ifdef PLATFORM_lapacke
lapack_int LAPACKE_zheev(int matrix_order, char jobz, char uplo, lapack_int n,
                         lapack_complex_double* a, lapack_int lda, double* w);
#elif defined PLATFORM_acml
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
void
daxpy_wrapper(LAPACK_INT n,        //number of elements of X,Y
              LAPACK_REAL alpha,   //scale factor
              const LAPACK_REAL* X, //pointer to head of vector X
              LAPACK_INT incx,     //increment with which to step through X
              LAPACK_REAL* Y,       //pointer to head of vector Y
              LAPACK_INT incy);     //increment with which to step through Y

//
// dnrm2
//
LAPACK_REAL
dnrm2_wrapper(LAPACK_INT N,
              const LAPACK_REAL* X,
              LAPACK_INT incx = 1);

//
// ddot
//
LAPACK_REAL
ddot_wrapper(LAPACK_INT N,
             const LAPACK_REAL* X,
             LAPACK_INT incx,
             const LAPACK_REAL* Y,
             LAPACK_INT incy);

//
// zdotc
//
Cplx
zdotc_wrapper(LAPACK_INT N,
              Cplx const* X,
              LAPACK_INT incx,
              Cplx const* Y,
              LAPACK_INT incy);

//
// dgemm
//
void
gemm_wrapper(bool transa, 
             bool transb,
             LAPACK_INT m,
             LAPACK_INT n,
             LAPACK_INT k,
             LAPACK_REAL alpha,
             LAPACK_REAL const* A,
             LAPACK_REAL const* B,
             LAPACK_REAL beta,
             LAPACK_REAL * C);

//
// zgemm
//
void
gemm_wrapper(bool transa, 
             bool transb,
             LAPACK_INT m,
             LAPACK_INT n,
             LAPACK_INT k,
             Cplx alpha,
             Cplx const* A,
             Cplx const* B,
             Cplx beta,
             Cplx * C);

//
// dgemv - matrix*vector multiply
//
void
gemv_wrapper(bool trans, 
             LAPACK_REAL alpha,
             LAPACK_REAL beta,
             LAPACK_INT m,
             LAPACK_INT n,
             const LAPACK_REAL* A,
             const LAPACK_REAL* x,
             LAPACK_INT incx,
             LAPACK_REAL* y,
             LAPACK_INT incy);

//
// zgemv - matrix*vector multiply
//
void
gemv_wrapper(bool trans, 
             Cplx alpha,
             Cplx beta,
             LAPACK_INT m,
             LAPACK_INT n,
             Cplx const* A,
             Cplx const* x,
             LAPACK_INT incx,
             Cplx* y,
             LAPACK_INT incy);


//
// dsyev
//
void
dsyev_wrapper(char jobz,        //if jobz=='V', compute eigs and evecs
              char uplo,        //if uplo=='U', read from upper triangle of A
              LAPACK_INT n,     //number of cols of A
              LAPACK_REAL* A,    //symmetric matrix A
              LAPACK_REAL* eigs, //eigenvalues on return
              LAPACK_INT& info);  //error info

//
// dscal
//
void
dscal_wrapper(LAPACK_INT N,
              LAPACK_REAL alpha,
              LAPACK_REAL* data,
              LAPACK_INT inc = 1);


void
dgesdd_wrapper(char * jobz,           //char* specifying how much of U, V to compute
                                    //choosing *jobz=='S' computes min(m,n) cols of U, V
               LAPACK_INT* m,       //number of rows of input matrix *A
               LAPACK_INT* n,       //number of cols of input matrix *A
               LAPACK_REAL *A,       //contents of input matrix A
               LAPACK_REAL *s,       //on return, singular values of A
               LAPACK_REAL *u,       //on return, unitary matrix U
               LAPACK_REAL *vt,      //on return, unitary matrix V transpose
               LAPACK_INT *info);

void
zgesdd_wrapper(char *jobz,           //char* specifying how much of U, V to compute
                                     //choosing *jobz=='S' computes min(m,n) cols of U, V
               LAPACK_INT *m,        //number of rows of input matrix *A
               LAPACK_INT *n,        //number of cols of input matrix *A
               Cplx *A,    //contents of input matrix A
               LAPACK_REAL *s,       //on return, singular values of A
               Cplx *u,    //on return, unitary matrix U
               Cplx *vt,   //on return, unitary matrix V transpose
               LAPACK_INT *info);


  void
dgesvd_wrapper(char * jobz,           //char* specifying how much of U, V to compute
                                    //choosing *jobz=='S' computes min(m,n) cols of U, V
               LAPACK_INT* m,       //number of rows of input matrix *A
               LAPACK_INT* n,       //number of cols of input matrix *A
               LAPACK_REAL *A,       //contents of input matrix A
               LAPACK_REAL *s,       //on return, singular values of A
               LAPACK_REAL *u,       //on return, unitary matrix U
               LAPACK_REAL *vt,      //on return, unitary matrix V transpose
               LAPACK_INT *info);

void
zgesvd_wrapper(char *jobz,           //char* specifying how much of U, V to compute
                                     //choosing *jobz=='S' computes min(m,n) cols of U, V
               LAPACK_INT *m,        //number of rows of input matrix *A
               LAPACK_INT *n,        //number of cols of input matrix *A
               Cplx *A,    //contents of input matrix A
               LAPACK_REAL *s,       //on return, singular values of A
               Cplx *u,    //on return, unitary matrix U
               Cplx *vt,   //on return, unitary matrix V transpose
               LAPACK_INT *info);


//
// dgeqrf
//
// QR factorization of a real matrix A
//
void
dgeqrf_wrapper(LAPACK_INT* m,     //number of rows of A
               LAPACK_INT* n,     //number of cols of A
               LAPACK_REAL* A,    //matrix A
                                  //on return upper triangle contains R
               LAPACK_INT* lda,   //size of A (usually same as n)
               LAPACK_REAL* tau,  //scalar factors of elementary reflectors
                                  //length should be min(m,n)
               LAPACK_INT* info);  //error info

//
// dorgqr
//
// Generates Q from output of QR factorization routine dgeqrf (see above)
//
void
dorgqr_wrapper(LAPACK_INT* m,     //number of rows of A
               LAPACK_INT* n,     //number of cols of A
               LAPACK_INT* k,     //number of elementary reflectors, typically min(m,n)
               LAPACK_REAL* A,    //matrix A, as returned from "A" argument of dgeqrf
                                  //on return contains Q
               LAPACK_INT* lda,   //size of A (usually same as n)
               LAPACK_REAL* tau,  //scalar factors as returned by dgeqrf
               LAPACK_INT* info);  //error info


  //
// dgeqrf
//
// QR factorization of a complex matrix A
//
void
zgeqrf_wrapper(LAPACK_INT* m,     //number of rows of A
               LAPACK_INT* n,     //number of cols of A
               Cplx* A,    //matrix A
                                  //on return upper triangle contains R
               LAPACK_INT* lda,   //size of A (usually same as n)
               LAPACK_COMPLEX* tau,  //scalar factors of elementary reflectors
                                  //length should be min(m,n)
               LAPACK_INT* info);  //error info

//
// dorgqr
//
// Generates Q from output of QR factorization routine zgeqrf (see above)
//
void
zungqr_wrapper(LAPACK_INT* m,     //number of rows of A
               LAPACK_INT* n,     //number of cols of A
               LAPACK_INT* k,     //number of elementary reflectors, typically min(m,n)
               Cplx* A,    //matrix A, as returned from "A" argument of dgeqrf
                                  //on return contains Q
               LAPACK_INT* lda,   //size of A (usually same as n)
               LAPACK_COMPLEX* tau,  //scalar factors as returned by zgeqrf
               LAPACK_INT* info);  //error info

// dgesv
//
// computes the solution to system of linear equations A*X = B
// where A is a general real matrix
//
LAPACK_INT
dgesv_wrapper(LAPACK_INT n,
              LAPACK_INT nrhs,
              LAPACK_REAL* a,
              LAPACK_REAL* b);

//
// zgesv
//
// computes the solution to system of linear euqations A*X =B
// where A is a general complex matrix
//
LAPACK_INT
zgesv_wrapper(LAPACK_INT n,
              LAPACK_INT nrhs,
              Cplx* a,
              Cplx* b);

//
// dlange
//
// returns the value of the 1-norm, Frobenius norm, infinity-norm, 
// or the largest absolute value of any element of a general rectangular matrix.
//
double
dlange_wrapper(char norm,
               LAPACK_INT m,
               LAPACK_INT n,
               double* a);

//
// zlange
//
// returns the value of the 1-norm, Frobenius norm, infinity-norm, 
// or the largest absolute value of any element of a general rectangular matrix.
//
LAPACK_REAL
zlange_wrapper(char norm,
               LAPACK_INT m,
               LAPACK_INT n,
               Cplx* a);

//
// zheev
//
// Eigenvalues and eigenvectors of complex Hermitian matrix A
//
LAPACK_INT 
zheev_wrapper(LAPACK_INT    N,  //number of cols of A
              Cplx        * A,  //matrix A, on return contains eigenvectors
              LAPACK_REAL * d); //eigenvalues on return

//
// dsygv
//
// Eigenvalues and eigenvectors of generalized eigenvalue problem
// A*x = lambda*B*x
// A and B must be symmetric
// B must be positive definite
//
void
dsygv_wrapper(char* jobz,           //if 'V', compute both eigs and evecs
                                    //if 'N', only eigenvalues
              char* uplo,           //if 'U', use upper triangle of A
              LAPACK_INT* n,        //number of cols of A
              LAPACK_REAL* A,       //matrix A, on return contains eigenvectors
              LAPACK_REAL* B,       //matrix B
              LAPACK_REAL* d,       //eigenvalues on return
              LAPACK_INT* info);  //error info

//
// dgeev
//
// Eigenvalues and eigenvectors of real, square matrix A
// A can be a general real matrix, not assumed symmetric
//
// Returns "info" integer
//
LAPACK_INT
dgeev_wrapper(char jobvl,          //if 'V', compute left eigenvectors, else 'N'
              char jobvr,          //if 'V', compute right eigenvectors, else 'N'
              LAPACK_INT n,        //number of rows/cols of A
              LAPACK_REAL const* A, //matrix A
              LAPACK_REAL* dr,      //real parts of eigenvalues
              LAPACK_REAL* di,      //imaginary parts of eigenvalues
              LAPACK_REAL* vl,      //left eigenvectors on return
              LAPACK_REAL* vr);     //right eigenvectors on return

//
// zgeev
//
// Eigenvalues and eigenvectors of complex, square matrix A
// A can be a general complex matrix, not assumed symmetric
//
// Returns "info" integer
//
LAPACK_INT
zgeev_wrapper(char jobvl,          //if 'V', compute left eigenvectors, else 'N'
              char jobvr,          //if 'V', compute right eigenvectors, else 'N'
              LAPACK_INT n,        //number of rows/cols of A
              Cplx const* A, //matrix A
              Cplx * d,    //eigenvalues
              Cplx * vl,   //left eigenvectors on return
              Cplx * vr);  //right eigenvectors on return

} //namespace itensor

#endif
