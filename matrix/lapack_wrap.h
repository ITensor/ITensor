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

#endif

//
// dsyev
//
int inline
dsyev_wrapper(char* jobz,
              char* uplo,
              int* n, 
              double* a,
              int* lda, 
              double* w, 
              double* work, 
              int* lwork,
              int* info)
    {
#ifdef PLATFORM_macos
    return dsyev_(jobz,uplo,n,a,lda,w,work,lwork,info);
#elif PLATFORM_acml
    return dsyev_(jobz,uplo,n,a,lda,w,work,lwork,info,1,1);
#endif
    }

#endif
