#ifndef __lapack_wrap_h
#define __lapack_wrap_h

//
// Headers and typedefs
//
#ifdef PLATFORM_macos

#include <Accelerate/Accelerate.h>
typedef __CLPK_integer
LAPACK_INT;
typedef __CLPK_doublereal
LAPACK_REAL;
typedef __CLPK_doublecomplex
LAPACK_COMPLEX;

#elif PLATFORM_acml

#include "acml.h"
typedef int
LAPACK_INT;
typedef double
LAPACK_REAL;
typedef doublecomplex
LAPACK_COMPLEX;

#elif PLATFORM_mkl

#include "mkl_lapack.h"
typedef MKL_INT
LAPACK_INT;
typedef double
LAPACK_REAL;
typedef MKL_Complex16
LAPACK_COMPLEX;

#endif


//
// dsyev
//
void inline
dsyev_wrapper(char* jobz,        //if jobz=='V', compute eigs and evecs
              char* uplo,        //if uplo=='U', read from upper triangle of A
              LAPACK_INT* n,     //number of cols of A
              LAPACK_REAL* A,    //symmetric matrix A
              LAPACK_INT* lda,   //size of A (usually same as n)
              LAPACK_REAL* eigs, //eigenvalues on return
              LAPACK_REAL* work, //workspace array
              LAPACK_INT* lwork, //size of workspace array
              LAPACK_INT* info)  //error info
    {
#ifdef PLATFORM_acml
    dsyev_(jobz,uplo,n,A,lda,eigs,work,lwork,info,1,1);
#else
    dsyev_(jobz,uplo,n,A,lda,eigs,work,lwork,info);
#endif
    }

void inline
zgesdd_wrapper(char *jobz,           //char* specifying how much of U, V to compute
                                     //choosing *jobz=='S' computes min(m,n) cols of U, V
               LAPACK_INT *m,        //number of rows of input matrix *A
               LAPACK_INT *n,        //number of cols of input matrix *A
               LAPACK_COMPLEX *A,    //contents of input matrix A
               LAPACK_INT *lda,      //typically lda==m
               LAPACK_REAL *s,       //on return, singular values of A
               LAPACK_COMPLEX *u,    //on return, unitary matrix U
               LAPACK_INT *ldu,      //typically ldu==m
               LAPACK_COMPLEX *vt,   //on return, unitary matrix V transpose
               LAPACK_INT *ldvt,     //typically ldvt==n
               LAPACK_COMPLEX *work, 
               LAPACK_INT *lwork, 
               LAPACK_REAL *rwork, 
               LAPACK_INT *iwork, 
               LAPACK_INT *info)
    {
    zgesdd_(jobz,m,n,A,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,iwork,info);
    }

#endif
