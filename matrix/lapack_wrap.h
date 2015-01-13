//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_lapack_wrap_h
#define __ITENSOR_lapack_wrap_h

//
// Headers and typedefs
//

#ifdef PLATFORM_lapack

#define LAPACK_REQUIRE_EXTERN

#include <complex>

namespace itensor {
typedef int
LAPACK_INT;
typedef double
LAPACK_REAL;
typedef struct
{
  double real, imag;
} LAPACK_COMPLEX;
};

#elif defined PLATFORM_macos

#include <Accelerate/Accelerate.h>
namespace itensor {
typedef __CLPK_integer
LAPACK_INT;
typedef __CLPK_doublereal
LAPACK_REAL;
typedef __CLPK_doublecomplex
LAPACK_COMPLEX;
};

#elif defined PLATFORM_acml

#define LAPACK_REQUIRE_EXTERN

//#include "acml.h"
namespace itensor {
typedef int
LAPACK_INT;
typedef double
LAPACK_REAL;
typedef struct
{
  double real, imag;
} LAPACK_COMPLEX;
};

#elif defined PLATFORM_mkl

#include "mkl_lapack.h"
namespace itensor {
typedef MKL_INT
LAPACK_INT;
typedef double
LAPACK_REAL;
typedef MKL_Complex16
LAPACK_COMPLEX;
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

#ifdef PLATFORM_acml
void F77NAME(dsyev)(char *jobz, char *uplo, int *n, double *a, int *lda, 
                    double *w, double *work, int *lwork, int *info, 
                    int jobz_len, int uplo_len);
#else
void F77NAME(dsyev)(const char* jobz, const char* uplo, const LAPACK_INT* n, double* a,
            const LAPACK_INT* lda, double* w, double* work, const LAPACK_INT* lwork,
            LAPACK_INT* info );
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

void F77NAME(dgehrd)(LAPACK_INT* n, LAPACK_INT* ilo, LAPACK_INT* ihi, LAPACK_REAL* A,
                     LAPACK_INT* lda, LAPACK_REAL* tau, LAPACK_INT* info);

void F77NAME(dorghr)(LAPACK_INT* n, LAPACK_INT* ilo, LAPACK_INT* ihi, LAPACK_REAL* A,
                     LAPACK_INT* lda, LAPACK_REAL* tau, LAPACK_INT* info);

void F77NAME(dhseqr)(LAPACK_INT* n, LAPACK_INT* ilo, LAPACK_INT* ihi, LAPACK_REAL* h,
                     LAPACK_INT* ldh, LAPACK_REAL* wr, LAPACK_REAL* wi, LAPACK_REAL* z,
                     LAPACK_INT* ldz, LAPACK_INT* info);

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
// dsyev
//
void inline
dsyev_wrapper(char* jobz,        //if jobz=='V', compute eigs and evecs
              char* uplo,        //if uplo=='U', read from upper triangle of A
              LAPACK_INT* n,     //numbec of cols of A
              LAPACK_REAL* A,    //symmetric matrix A
              LAPACK_INT* lda,   //size of A (usually same as n)
              LAPACK_REAL* eigs, //eigenvalues on return
              LAPACK_INT* info)  //error info
    {
    LAPACK_INT lwork = max(1,3*(*n)-1);
    LAPACK_REAL work[lwork];

#ifdef PLATFORM_acml
    F77NAME(dsyev)(jobz,uplo,n,A,lda,eigs,work,&lwork,info,1,1);
#else
    F77NAME(dsyev)(jobz,uplo,n,A,lda,eigs,work,&lwork,info);
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
    LAPACK_INT l = min(*m,*n),
               g = max(*m,*n);
    LAPACK_INT lwork = l*l+2*l+g+100;
    LAPACK_COMPLEX work[lwork];
    LAPACK_REAL rwork[5*l*(1+l)];
    LAPACK_INT iwork[8*l];
#ifdef PLATFORM_acml
    LAPACK_INT jobz_len = 1;
    F77NAME(zgesdd)(jobz,m,n,A,m,s,u,m,vt,n,work,&lwork,rwork,iwork,info,jobz_len);
#else
    F77NAME(zgesdd)(jobz,m,n,A,m,s,u,m,vt,n,work,&lwork,rwork,iwork,info);
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
    int lwork = max(1,4*max(*n,*m));
    LAPACK_REAL work[lwork]; 
    F77NAME(dgeqrf)(m,n,A,lda,tau,work,&lwork,info);
    }

//
// dgehrd
//
// Reduces a general matrix to upper Hessenberg form.
//
void inline
dgehrd_wrapper(LAPACK_INT* n,     //order of the matrix A
               LAPACK_INT* ilo,   //
               LAPACK_INT* ihi,   //it is assumed that A is already upper triangular in rows 
                                  //and columns 1:ilo-1 and ihi+1:n 
                                  //ilo and ihi are normally set by a previous call to dgebal; 
                                  //otherwise they should be set to 1 and n respectively
               LAPACK_REAL* A,    //on entry, the n-by-n general matrix to be reduced
                                  //on exit, the upper triangle and the first subdiagonal of A
                                  //are overwritten with the upper Hessenberg matrix H, and the
                                  //elements below the first subdiagonal, with the array tau,
                                  //represent the orthogonal matrix Q as a product of elementary
                                  //reflectors
               LAPACK_INT* lda,   //The leading dimension of the array A
                                  //lda >= max(1,n)
               LAPACK_REAL* tau,  //the scalar factors of the elementary reflectors
                                  //elements 1:ilo-1 and ihi:n-1 of tau are set to zero.
               LAPACK_INT* info)  //error info
    {
    int lwork = max(1,4*max(*n,*n));
    LAPACK_REAL work[lwork];
    F77NAME(dgehrd)(n,ilo,ihi,A,lda,tau,work,&lwork,info);
    }

//
// dorghr
//
// Generates the real orthogonal matrix Q determined by dgehrd.
//
void inline
dorghr_wrapper(LAPACK_INT* n,     //order of the matrix Q
               LAPACK_INT* ilo,   //
               LAPACK_INT* ihi,   //ilo and ihi must have the same values as in the previous call
                                  //of dgehrd. Q is equal to the unit matrix except in the
                                  //submatrix Q(ilo+1:ihi,ilo+1:ihi).
               LAPACK_REAL* A,    //on entry, the vectors which define the elementary reflectors,
                                  //as returned by dgehrd
                                  //on exit, the n-by-n orthogonal matrix Q
               LAPACK_INT* lda,   //the leading dimension of the array A
                                  //lda >= max(1,n)
               LAPACK_REAL* tau,  //tau(i) must contain the scalar factor of the elementary
                                  //reflector H(i), as returned by dgehrd
               LAPACK_INT* info)  //error info
    {
    int lwork = max(1,4*max(*n,*n));
    LAPACK_REAL work[lwork];
    F77NAME(dorghr)(n,ilo,ihi,A,lda,tau,work,&lwork,info);
    }

//
// dhseqr
//
// Computes all eigenvalues and (optionally) the Schur factorization 
// of a matrix reduced to Hessenberg form.
//
void inline
dhseqr_wrapper(LAPACK_INT* n,     //order of the matrix H
               LAPACK_INT* ilo,   //
               LAPACK_INT* ihi,   //it is assumed that H is already upper triangular in rows
                                  //and columns 1:ilo-1 and ihi+1:n 
                                  //ilo and ihi are normally set by a previous call to dgebal, 
                                  //and then passed to zgehrd when the matrix output by dgebal 
                                  //is reduced to Hessenberg form
                                  //otherwise ilo and ihi should be set to 1 and n respectively
               LAPACK_REAL* h,    //on entry, the upper Hessenberg matrix H
                                  //on exit, if info = 0 and job = 'S', then H contains the
                                  //upper quasi-triangular matrix T from the Schur decomposition
                                  //(the Schur form); 2-by-2 diagonal blocks (corresponding to
                                  //complex conjugate pairs of eigenvalues) are returned in
                                  //standard form, with H(i,i) = H(i+1,i+1) and
                                  //H(i+1,i)*H(i,i+1) <= 0
                                  //if info = 0 and job = 'E', the contents of H are unspecified 
                                  //on exit
               LAPACK_INT* ldh,   //leading dimension of the array H
                                  //ldh >= max(1,n)
               LAPACK_REAL* wr,   //
               LAPACK_REAL* wi,   //the real and imaginary parts, respectively, of the computed
                                  //eigenvalues 
                                  //if two eigenvalues are computed as a complex conjugate pair, 
                                  //they are stored in consecutive elements of wr and wi, say the 
                                  //i-th and (i+1)th, with wi(i) > 0 and wi(i+1) < 0 
                                  //if job = 'S', the eigenvalues are stored in the same order as 
                                  //on the diagonal of the Schur form returned in H, with 
                                  //wr(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2 diagonal block, 
                                  //wi(i) = sqrt(-H(i+1,i)*H(i,i+1)) and wi(i+1) = -wi(i)
               LAPACK_REAL* z,    //if compz = 'N', Z is not referenced
                                  //if compz = 'I', on entry Z need not be set and on exit,
                                  //if info = 0, Z contains the orthogonal matrix Z of the Schur
                                  //vectors of H
                                  //if compz = 'V', on entry Z must contain an
                                  //n-by-n matrix Q, which is assumed to be equal to the unit
                                  //matrix except for the submatrix Z(ILO:IHI,ILO:IHI)
                                  //on exit, if info = 0, Z contains Q*Z
                                  //normally Q is the orthogonal matrix generated by dorghr
                                  //after the call to dgehrd which formed the Hessenberg matrix H
               LAPACK_INT* ldz,   //the leading dimension of the array Z
                                  //if compz = 'I' or compz = 'V', then ldz >= max(1,n)  
                                  //otherwize, ldz >= 1
               LAPACK_INT* info)  //error info
    {
    int lwork = max(1,4*max(*n,*n));
    LAPACK_REAL work[lwork];
    char job = 'S';
    char compz = 'V';
    F77NAME(dhseqr)(&job,&compz,n,ilo,ihi,h,ldh,wr,wi,z,ldz,work,&lwork,info);
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
    int lwork = max(1,4*max(*n,*m));
    LAPACK_REAL work[lwork]; 
    F77NAME(dorgqr)(m,n,k,A,lda,tau,work,&lwork,info);
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
    int itype = 1;
    LAPACK_INT lwork = max(1,3*(*n)-1);//max(1, 1+6*N+2*N*N);
    LAPACK_REAL work[lwork];
#ifdef PLATFORM_acml
    LAPACK_INT jobz_len = 1;
    LAPACK_INT uplo_len = 1;
    F77NAME(dsygv)(&itype,jobz,uplo,n,A,n,B,n,d,work,&lwork,info,jobz_len,uplo_len);
#else
    F77NAME(dsygv)(&itype,jobz,uplo,n,A,n,B,n,d,work,&lwork,info);
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
    int nevecl = (*jobvl == 'V' ? *n : 1);
    int nevecr = (*jobvr == 'V' ? *n : 1);
    LAPACK_INT lwork = max(1,4*(*n));
    LAPACK_REAL work[lwork];
#ifdef PLATFORM_acml
    LAPACK_INT jobvl_len = 1;
    LAPACK_INT jobvr_len = 1;
    F77NAME(dgeev)(jobvl,jobvr,n,A,n,dr,di,vl,&nevecl,vr,&nevecr,work,&lwork,info,jobvl_len,jobvr_len);
#else
    F77NAME(dgeev)(jobvl,jobvr,n,A,n,dr,di,vl,&nevecl,vr,&nevecr,work,&lwork,info);
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
    int nevecl = (*jobvl == 'V' ? *n : 1);
    int nevecr = (*jobvr == 'V' ? *n : 1);
    LAPACK_INT lwork = max(1,4*(*n));
    LAPACK_COMPLEX work[lwork];
    LAPACK_INT lrwork = max(1,2*(*n));
    LAPACK_REAL rwork[lrwork];
#ifdef PLATFORM_acml
    F77NAME(zgeev)(jobvl,jobvr,n,A,n,d,vl,&nevecl,vr,&nevecr,work,&lwork,rwork,info,1,1);
#else
    F77NAME(zgeev)(jobvl,jobvr,n,A,n,d,vl,&nevecl,vr,&nevecr,work,&lwork,rwork,info);
#endif
    }

}; //namespace itensor

#endif
