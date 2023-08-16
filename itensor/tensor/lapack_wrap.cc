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
#include "itensor/tensor/lapack_wrap.h"
//#include "itensor/tensor/permutecplx.h"

namespace itensor {

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
              LAPACK_INT incy)     //increment with which to step through Y
    {
#ifdef ITENSOR_USE_CBLAS
    cblas_daxpy(n,alpha,X,incx,Y,incy);
#else
    auto Xnc = const_cast<LAPACK_REAL*>(X);
    F77NAME(daxpy)(&n,&alpha,Xnc,&incx,Y,&incy);
#endif
    }

//
// dnrm2
//
LAPACK_REAL 
dnrm2_wrapper(LAPACK_INT N,
              const LAPACK_REAL* X,
              LAPACK_INT incx)
    {
#ifdef ITENSOR_USE_CBLAS
    return cblas_dnrm2(N,X,incx);
#else
    auto *Xnc = const_cast<LAPACK_REAL*>(X);
    return F77NAME(dnrm2)(&N,Xnc,&incx);
#endif
    return -1;
    }

//
// ddot
//
LAPACK_REAL 
ddot_wrapper(LAPACK_INT N,
             const LAPACK_REAL* X,
             LAPACK_INT incx,
             const LAPACK_REAL* Y,
             LAPACK_INT incy)
    {
#ifdef ITENSOR_USE_CBLAS
    return cblas_ddot(N,X,incx,Y,incy);
#else
    auto *Xnc = const_cast<LAPACK_REAL*>(X);
    auto *Ync = const_cast<LAPACK_REAL*>(Y);
    return F77NAME(ddot)(&N,Xnc,&incx,Ync,&incy);
#endif
    return -1;
    }

//
// zdotc
//
Cplx
zdotc_wrapper(LAPACK_INT N,
              Cplx const* X,
              LAPACK_INT incx,
              Cplx const* Y,
              LAPACK_INT incy)
    {
#ifdef ITENSOR_USE_CBLAS
    Cplx res;
#if defined PLATFORM_openblas
    auto pX = reinterpret_cast<OPENBLAS_CONST double*>(X);
    auto pY = reinterpret_cast<OPENBLAS_CONST double*>(Y);
    auto pres = reinterpret_cast<openblas_complex_double*>(&res);
#else
    auto pX = reinterpret_cast<const void*>(X);
    auto pY = reinterpret_cast<const void*>(Y);
    auto pres = reinterpret_cast<void*>(&res);
#endif
    cblas_zdotc_sub(N,pX,incx,pY,incy,pres);
    return res;
#else
    auto ncX = const_cast<Cplx*>(X);
    auto ncY = const_cast<Cplx*>(Y);
    auto pX = reinterpret_cast<LAPACK_COMPLEX*>(ncX);
    auto pY = reinterpret_cast<LAPACK_COMPLEX*>(ncY);
    auto res = F77NAME(zdotc)(&N,pX,&incx,pY,&incy);
    auto cplx_res = reinterpret_cast<Cplx*>(&res);
    return *cplx_res;
#endif
    return Cplx{};
    }

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
             LAPACK_REAL * C)
    {
    LAPACK_INT lda = m,
               ldb = k;
#ifdef ITENSOR_USE_CBLAS
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
// zgemm
//
void 
gemm_wrapper(bool transa, 
             bool transb,
             LAPACK_INT m,
             LAPACK_INT n,
             LAPACK_INT k,
             Cplx alpha,
             const Cplx* A,
             const Cplx* B,
             Cplx beta,
             Cplx* C)
    {
    LAPACK_INT lda = m,
               ldb = k;
#ifdef PLATFORM_openblas
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
    //auto ralpha = realRef(alpha);
    //auto ialpha = imagRef(alpha);
    //auto rbeta = realRef(beta);
    //auto ibeta = imagRef(beta);
    //if(ialpha != 0.0 || ibeta != 0.0)
    //    {
    //    throw std::runtime_error("Complex alpha, beta not supported in zgemm for PLATFORM=openblas");
    //    }
    auto* palpha = reinterpret_cast<double*>(&alpha);
    auto* pbeta = reinterpret_cast<double*>(&beta);
    auto* pA = reinterpret_cast<const double*>(A);
    auto* pB = reinterpret_cast<const double*>(B);
    auto* pC = reinterpret_cast<double*>(C);
	cblas_zgemm(CblasColMajor,at,bt,m,n,k,palpha,pA,lda,pB,ldb,pbeta,pC,m);
#else //platform not openblas
#ifdef ITENSOR_USE_CBLAS
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
    auto palpha = (void*)(&alpha); 
    auto pbeta = (void*)(&beta); 
    cblas_zgemm(CblasColMajor,at,bt,m,n,k,palpha,(void*)A,lda,(void*)B,ldb,pbeta,(void*)C,m);
#else //use Fortran zgemm
    auto *ncA = const_cast<Cplx*>(A);
    auto *ncB = const_cast<Cplx*>(B);
    auto *pA = reinterpret_cast<LAPACK_COMPLEX*>(ncA);
    auto *pB = reinterpret_cast<LAPACK_COMPLEX*>(ncB);
    auto *pC = reinterpret_cast<LAPACK_COMPLEX*>(C);
    auto *palpha = reinterpret_cast<LAPACK_COMPLEX*>(&alpha);
    auto *pbeta = reinterpret_cast<LAPACK_COMPLEX*>(&beta);
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
    F77NAME(zgemm)(&at,&bt,&m,&n,&k,palpha,pA,&lda,pB,&ldb,pbeta,pC,&m);
#endif
#endif
    }

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
             LAPACK_INT incy)
    {
#ifdef ITENSOR_USE_CBLAS
    auto Tr = trans ? CblasTrans : CblasNoTrans;
    cblas_dgemv(CblasColMajor,Tr,m,n,alpha,A,m,x,incx,beta,y,incy);
#else
    char Tr = trans ? 'T' : 'N';
    F77NAME(dgemv)(&Tr,&m,&n,&alpha,const_cast<LAPACK_REAL*>(A),&m,const_cast<LAPACK_REAL*>(x),&incx,&beta,y,&incy);
#endif
    }

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
             LAPACK_INT incy)
    {
#ifdef PLATFORM_openblas
    auto Tr = trans ? CblasTrans : CblasNoTrans;
    //auto ralpha = realRef(alpha);
    //auto ialpha = imagRef(alpha);
    //auto rbeta = realRef(beta);
    //auto ibeta = imagRef(beta);
    //if(ialpha != 0.0 || ibeta != 0.0)
    //    {
    //    throw std::runtime_error("Complex alpha, beta not supported in zgemm for PLATFORM=openblas");
    //    }
	auto* palpha = reinterpret_cast<double*>(&alpha);
	auto* pbeta = reinterpret_cast<double*>(&beta);
	auto* pA = reinterpret_cast<const double*>(A);
	auto* px = reinterpret_cast<const double*>(x);
	auto* py = reinterpret_cast<double*>(y);
    cblas_zgemv(CblasColMajor,Tr,m,n,palpha,pA,m,px,incx,pbeta,py,incy);
#else //platform other than openblas
#ifdef ITENSOR_USE_CBLAS
    auto Tr = trans ? CblasTrans : CblasNoTrans;
    auto palpha = reinterpret_cast<void*>(&alpha); 
    auto pbeta = reinterpret_cast<void*>(&beta); 
    cblas_zgemv(CblasColMajor,Tr,m,n,palpha,(void*)A,m,(void*)x,incx,pbeta,(void*)y,incy);
#else
    char Tr = trans ? 'T' : 'N';
    auto ncA = const_cast<Cplx*>(A); 
    auto ncx = const_cast<Cplx*>(x); 
    auto pA = reinterpret_cast<LAPACK_COMPLEX*>(ncA); 
    auto px = reinterpret_cast<LAPACK_COMPLEX*>(ncx); 
    auto py = reinterpret_cast<LAPACK_COMPLEX*>(y); 
    auto palpha = reinterpret_cast<LAPACK_COMPLEX*>(&alpha); 
    auto pbeta = reinterpret_cast<LAPACK_COMPLEX*>(&beta); 
    F77NAME(zgemv)(&Tr,&m,&n,palpha,pA,&m,px,&incx,pbeta,py,&incy);
#endif
#endif
    }


//
// dsyev
//
void 
dsyev_wrapper(char jobz,        //if jobz=='V', compute eigs and evecs
              char uplo,        //if uplo=='U', read from upper triangle of A
              LAPACK_INT n,     //number of cols of A
              LAPACK_REAL* A,    //symmetric matrix A
              LAPACK_REAL* eigs, //eigenvalues on return
              LAPACK_INT& info)  //error info
    {
    std::vector<LAPACK_REAL> work;
    LAPACK_INT lda = n;

#ifdef PLATFORM_acml
    static const LAPACK_INT one = 1;
    LAPACK_INT lwork = std::max(one,3*n-1);
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
void 
dscal_wrapper(LAPACK_INT N,
              LAPACK_REAL alpha,
              LAPACK_REAL* data,
              LAPACK_INT inc)
    {
#ifdef ITENSOR_USE_CBLAS
    cblas_dscal(N,alpha,data,inc);
#else
    F77NAME(dscal)(&N,&alpha,data,&inc);
#endif
    }

void 
zgesdd_wrapper(char *jobz,           //char* specifying how much of U, V to compute
                                     //choosing *jobz=='S' computes min(m,n) cols of U, V
               LAPACK_INT *m,        //number of rows of input matrix *A
               LAPACK_INT *n,        //number of cols of input matrix *A
               Cplx *A,    //contents of input matrix A
               LAPACK_REAL *s,       //on return, singular values of A
               Cplx *u,    //on return, unitary matrix U
               Cplx *vt,   //on return, unitary matrix V transpose
               LAPACK_INT *info)
    {
    std::vector<LAPACK_COMPLEX> work;
    std::vector<LAPACK_REAL> rwork;
    std::vector<LAPACK_INT> iwork;
    auto pA = reinterpret_cast<LAPACK_COMPLEX*>(A);
    auto pU = reinterpret_cast<LAPACK_COMPLEX*>(u);
    auto pVt = reinterpret_cast<LAPACK_COMPLEX*>(vt);
    LAPACK_INT l = std::min(*m,*n),
               g = std::max(*m,*n);
    LAPACK_INT lwork = l*l+2*l+g+100;
    work.resize(lwork);
    rwork.resize(5*l*(1+l));
    iwork.resize(8*l);
#ifdef PLATFORM_acml
    LAPACK_INT jobz_len = 1;
    F77NAME(zgesdd)(jobz,m,n,pA,m,s,pU,m,pVt,&l,work.data(),&lwork,rwork.data(),iwork.data(),info,jobz_len);
#else
    F77NAME(zgesdd)(jobz,m,n,pA,m,s,pU,m,pVt,&l,work.data(),&lwork,rwork.data(),iwork.data(),info);
#endif
    }



  void 
dgesdd_wrapper(char* jobz,           //char* specifying how much of U, V to compute
                                    //choosing *jobz=='S' computes min(m,n) cols of U, V
               LAPACK_INT *m,        //number of rows of input matrix *A
               LAPACK_INT *n,        //number of cols of input matrix *A
               LAPACK_REAL *A,      //contents of input matrix A
               LAPACK_REAL *s,      //on return, singular values of A
               LAPACK_REAL *u,           //on return, unitary matrix U
               LAPACK_REAL *vt,          //on return, unitary matrix V transpose
               LAPACK_INT *info)
    {
    std::vector<LAPACK_REAL> work;
    std::vector<LAPACK_INT> iwork;
    LAPACK_INT l = std::min(*m,*n),
               g = std::max(*m,*n);
    LAPACK_INT lwork = l*(6 + 4*l) + g;
    work.resize(lwork);
    iwork.resize(8*l);
#ifdef PLATFORM_acml
    LAPACK_INT jobz_len = 1;
    F77NAME(dgesdd)(jobz,m,n,A,m,s,u,m,vt,&l,work.data(),&lwork,iwork.data(),info,jobz_len);
#else
    F77NAME(dgesdd)(jobz,m,n,A,m,s,u,m,vt,&l,work.data(),&lwork,iwork.data(),info);
#endif
    }



  void 
zgesvd_wrapper(char *jobz,           //char* specifying how much of U, V to compute
                                     //choosing *jobz=='S' computes min(m,n) cols of U, V
               LAPACK_INT *m,        //number of rows of input matrix *A
               LAPACK_INT *n,        //number of cols of input matrix *A
               Cplx *A,    //contents of input matrix A
               LAPACK_REAL *s,       //on return, singular values of A
               Cplx *u,    //on return, unitary matrix U
               Cplx *vt,   //on return, unitary matrix V transpose
               LAPACK_INT *info)
    {
    std::vector<LAPACK_COMPLEX> work;
    std::vector<LAPACK_REAL> rwork;
    std::vector<LAPACK_INT> iwork;
    auto pA = reinterpret_cast<LAPACK_COMPLEX*>(A);
    auto pU = reinterpret_cast<LAPACK_COMPLEX*>(u);
    auto pVt = reinterpret_cast<LAPACK_COMPLEX*>(vt);
    LAPACK_INT l = std::min(*m,*n),
               g = std::max(*m,*n);
    LAPACK_INT lwork = l*l+2*l+g+100;
    work.resize(lwork);
    rwork.resize(5*l*(1+l));
    iwork.resize(8*l);
#ifdef PLATFORM_acml
    LAPACK_INT jobz_len = 1;
    F77NAME(zgesvd)(jobz,jobz,m,n,pA,m,s,pU,m,pVt,&l,work.data(),&lwork,rwork.data(),info,jobz_len);
#else
    F77NAME(zgesvd)(jobz,jobz,m,n,pA,m,s,pU,m,pVt,&l,work.data(),&lwork,rwork.data(),info);
#endif
    }



  void 
dgesvd_wrapper(char* jobz,           //char* specifying how much of U, V to compute
                                    //choosing *jobz=='S' computes min(m,n) cols of U, V
               LAPACK_INT *m,        //number of rows of input matrix *A
               LAPACK_INT *n,        //number of cols of input matrix *A
               LAPACK_REAL *A,      //contents of input matrix A
               LAPACK_REAL *s,      //on return, singular values of A
               LAPACK_REAL *u,           //on return, unitary matrix U
               LAPACK_REAL *vt,          //on return, unitary matrix V transpose
               LAPACK_INT *info)
    {
    std::vector<LAPACK_REAL> work;
    // std::vector<LAPACK_REAL> superb;
    
    std::vector<LAPACK_INT> iwork;
    LAPACK_INT l = std::min(*m,*n),
               g = std::max(*m,*n);
    LAPACK_INT lwork = l*(6 + 4*l) + g;
    work.resize(lwork);
    iwork.resize(8*l);
    //superb.resize(l -1);
#ifdef PLATFORM_acml
    LAPACK_INT jobz_len = 1;
    F77NAME(dgesvd)(jobz,jobz,m,n,A,m,s,u,m,vt,&l,work.data(),&lwork, info, jobz_len);
#else
    F77NAME(dgesvd)(jobz,jobz,m,n,A,m,s,u,m,vt,&l,work.data(),&lwork, info);
#endif
    }

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
               LAPACK_INT* info)  //error info
    {
    static const LAPACK_INT one = 1;
    std::vector<LAPACK_REAL> work;
    LAPACK_INT lwork = std::max(one,4*std::max(*n,*m));
    work.resize(lwork+2); 
    F77NAME(dgeqrf)(m,n,A,lda,tau,work.data(),&lwork,info);
    }

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
               LAPACK_INT* info)  //error info
    {
    static const LAPACK_INT one = 1;
    std::vector<LAPACK_REAL> work;
    auto lwork = std::max(one,4*std::max(*n,*m));
    work.resize(lwork+2); 
    F77NAME(dorgqr)(m,n,k,A,lda,tau,work.data(),&lwork,info);
    }


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
               LAPACK_INT* info)  //error info
    {
    static const LAPACK_INT one = 1;
    std::vector<LAPACK_COMPLEX> work;
    LAPACK_INT lwork = std::max(one,4*std::max(*n,*m));
    work.resize(lwork+2);
    static_assert(sizeof(LAPACK_COMPLEX)==sizeof(Cplx),"LAPACK_COMPLEX and itensor::Cplx have different size");
    auto pA = reinterpret_cast<LAPACK_COMPLEX*>(A);
    F77NAME(zgeqrf)(m,n,pA,lda,tau,work.data(),&lwork,info);
    }

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
               LAPACK_COMPLEX* tau,  //scalar factors as returned by dgeqrf
               LAPACK_INT* info)  //error info
    {
    static const LAPACK_INT one = 1;
    std::vector<LAPACK_COMPLEX> work;
    auto lwork = std::max(one,4*std::max(*n,*m));
    work.resize(lwork+2);
    static_assert(sizeof(LAPACK_COMPLEX)==sizeof(Cplx),"LAPACK_COMPLEX and itensor::Cplx have different size");
    auto pA = reinterpret_cast<LAPACK_COMPLEX*>(A);
    #ifdef PLATFORM_lapacke
    LAPACKE_zungqr(LAPACK_COL_MAJOR,jobz,uplo,N,A,N,w.data());
    #else
    F77NAME(zungqr)(m,n,k,pA,lda,tau,work.data(),&lwork,info);
    #endif
    }

//
// dgesv
//
LAPACK_INT
dgesv_wrapper(LAPACK_INT n,
              LAPACK_INT nrhs,
              LAPACK_REAL* a,
              LAPACK_REAL* b)
	{
	LAPACK_INT lda = n;
	std::vector<LAPACK_INT> ipiv(n);
	LAPACK_INT ldb = n;
	LAPACK_INT info = 0;
	F77NAME(dgesv)(&n,&nrhs,a,&lda,ipiv.data(),b,&ldb,&info);
	return info;
	}

//
// zgesv
//
LAPACK_INT
zgesv_wrapper(LAPACK_INT n,
              LAPACK_INT nrhs,
              Cplx* a,
              Cplx* b)
	{
	auto pa = reinterpret_cast<LAPACK_COMPLEX*>(a);
	auto pb = reinterpret_cast<LAPACK_COMPLEX*>(b);
	LAPACK_INT lda = n;
	std::vector<LAPACK_INT> ipiv(n);
	LAPACK_INT ldb = n;
	LAPACK_INT info = 0;
	F77NAME(zgesv)(&n,&nrhs,pa,&lda,ipiv.data(),pb,&ldb,&info);
	return info;
	}

//
// dlange
//
double
dlange_wrapper(char norm,
               LAPACK_INT m,
               LAPACK_INT n,
               double* a)
	{
	double norma;
#ifdef PLATFORM_lapacke
	norma = LAPACKE_dlange(LAPACK_COL_MAJOR,norm,m,n,a,m);
#else
	std::vector<double> work;
	if(norm == 'I' || norm == 'i') work.resize(m);
#ifdef PLATFORM_acml
	LAPACK_INT norm_len = 1;
	norma = F77NAME(dlange)(&norm,&m,&n,a,&m,work.data(),norm_len);
#else
	norma = F77NAME(dlange)(&norm,&m,&n,a,&m,work.data());
#endif
#endif
	return norma;
	}

//
// zlange
//
LAPACK_REAL
zlange_wrapper(char norm,
               LAPACK_INT m,
               LAPACK_INT n,
               Cplx* a)
	{
	LAPACK_REAL norma;
#ifdef PLATFORM_lapacke
	auto pA = reinterpret_cast<lapack_complex_double*>(a);
	norma = LAPACKE_zlange(LAPACK_COL_MAJOR,norm,m,n,pa,m);
#else
	std::vector<double> work;
	if(norm == 'I' || norm == 'i') work.resize(m);
	auto pA = reinterpret_cast<LAPACK_COMPLEX*>(a);
#ifdef PLATFORM_acml
	LAPACK_INT norm_len = 1;
	norma = F77NAME(zlange)(&norm,&m,&n,pA,&m,work.data(),norm_len);
#else
	norma = F77NAME(zlange)(&norm,&m,&n,pA,&m,work.data());
#endif
#endif
	return norma;
	}

//
// zheev
//
// Eigenvalues and eigenvectors of complex Hermitian matrix A
//
LAPACK_INT 
zheev_wrapper(LAPACK_INT      N,  //number of cols of A
              Cplx          * A,  //matrix A, on return contains eigenvectors
              LAPACK_REAL   * d)  //eigenvalues on return
    {
    static const LAPACK_INT one = 1;
    char jobz = 'V';
    char uplo = 'U';
#ifdef PLATFORM_lapacke
    std::vector<LAPACK_REAL> work(N);
    LAPACKE_zheev(LAPACK_COL_MAJOR,jobz,uplo,N,A,N,w.data());
#else
    LAPACK_INT lwork = std::max(one,3*N-1);//max(1, 1+6*N+2*N*N);
    std::vector<LAPACK_COMPLEX> work(lwork);
    std::vector<LAPACK_REAL> rwork(lwork);
    LAPACK_INT info = 0;
    static_assert(sizeof(LAPACK_COMPLEX)==sizeof(Cplx),"LAPACK_COMPLEX and itensor::Cplx have different size");
    auto pA = reinterpret_cast<LAPACK_COMPLEX*>(A);
#ifdef PLATFORM_acml
    LAPACK_INT jobz_len = 1;
    LAPACK_INT uplo_len = 1;
    F77NAME(zheev)(&jobz,&uplo,&N,pA,&N,d,work.data(),&lwork,rwork.data(),&info,jobz_len,uplo_len);
#else
    F77NAME(zheev)(&jobz,&uplo,&N,pA,&N,d,work.data(),&lwork,rwork.data(),&info);
#endif

#endif //PLATFORM_lapacke
    return info;
    }

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
              LAPACK_INT* info)  //error info
    {
    static const LAPACK_INT one = 1;
    std::vector<LAPACK_REAL> work;
    LAPACK_INT itype = 1;
    LAPACK_INT lwork = std::max(one,3*(*n)-1);//std::max(1, 1+6*N+2*N*N);
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
LAPACK_INT 
dgeev_wrapper(char jobvl,          //if 'V', compute left eigenvectors, else 'N'
              char jobvr,          //if 'V', compute right eigenvectors, else 'N'
              LAPACK_INT n,        //number of rows/cols of A
              LAPACK_REAL const* A, //matrix A
              LAPACK_REAL* dr,      //real parts of eigenvalues
              LAPACK_REAL* di,      //imaginary parts of eigenvalues
              LAPACK_REAL* vl,      //left eigenvectors on return
              LAPACK_REAL* vr)      //right eigenvectors on return
    {
    std::vector<LAPACK_REAL> work;
    std::vector<LAPACK_REAL> cpA;

    cpA.resize(n*n);
    std::copy(A,A+n*n,cpA.data());
    
    LAPACK_INT nevecl = (jobvl == 'V' ? n : 1);
    LAPACK_INT nevecr = (jobvr == 'V' ? n : 1);
    LAPACK_INT info = 0;
#ifdef PLATFORM_acml
    LAPACK_INT lwork = -1;
    LAPACK_REAL wquery = 0;
    F77NAME(dgeev)(&jobvl,&jobvr,&n,cpA.data(),&n,dr,di,vl,&nevecl,vr,&nevecr,&wquery,&lwork,&info,1,1);

    lwork = static_cast<LAPACK_INT>(wquery);
    work.resize(lwork);
    F77NAME(dgeev)(&jobvl,&jobvr,&n,cpA.data(),&n,dr,di,vl,&nevecl,vr,&nevecr,work.data(),&lwork,&info,1,1);
#else
    LAPACK_INT lwork = -1;
    LAPACK_REAL wquery = 0;
    F77NAME(dgeev)(&jobvl,&jobvr,&n,cpA.data(),&n,dr,di,vl,&nevecl,vr,&nevecr,&wquery,&lwork,&info);

    lwork = static_cast<LAPACK_INT>(wquery);
    work.resize(lwork);
    F77NAME(dgeev)(&jobvl,&jobvr,&n,cpA.data(),&n,dr,di,vl,&nevecl,vr,&nevecr,work.data(),&lwork,&info);
#endif
    //println("jobvl = ",jobvl);
    //println("nevecl = ",nevecl);
    //println("vl data = ");
    //for(auto j = 0; j < n*n; ++j)
    //    {
    //    println(*vl);
    //    ++vl;
    //    }
    //println("vr data = ");
    //for(auto j = 0; j < n*n; ++j)
    //    {
    //    println(*vr);
    //    ++vr;
    //    }
    return info;
    }

//
// zgeev
//
// Eigenvalues and eigenvectors of complex, square matrix A
// A can be a general complex matrix, not assumed symmetric
//
LAPACK_INT 
zgeev_wrapper(char jobvl,          //if 'V', compute left eigenvectors, else 'N'
              char jobvr,          //if 'V', compute right eigenvectors, else 'N'
              LAPACK_INT n,        //number of rows/cols of A
              Cplx const* A, //matrix A
              Cplx * d,    //eigenvalues
              Cplx * vl,   //left eigenvectors on return
              Cplx * vr)   //right eigenvectors on return
    {
    static const LAPACK_INT one = 1;
    std::vector<LAPACK_COMPLEX> cpA;
    std::vector<LAPACK_COMPLEX> work;
    std::vector<LAPACK_REAL> rwork;
    LAPACK_INT nevecl = (jobvl == 'V' ? n : 1);
    LAPACK_INT nevecr = (jobvr == 'V' ? n : 1);
    LAPACK_INT lwork = std::max(one,4*n);
    work.resize(lwork);
    LAPACK_INT lrwork = std::max(one,2*n);
    rwork.resize(lrwork);

    //Copy A data into cpA
    cpA.resize(n*n);
    auto pA = reinterpret_cast<LAPACK_COMPLEX const*>(A);
    std::copy(pA,pA+n*n,cpA.data());

    auto pd = reinterpret_cast<LAPACK_COMPLEX*>(d);
    auto pvl = reinterpret_cast<LAPACK_COMPLEX*>(vl);
    auto pvr = reinterpret_cast<LAPACK_COMPLEX*>(vr);

    LAPACK_INT info = 0;
#ifdef PLATFORM_acml
    F77NAME(zgeev)(&jobvl,&jobvr,&n,cpA.data(),&n,pd,pvl,&nevecl,pvr,&nevecr,work.data(),&lwork,rwork.data(),&info,1,1);
#else
    F77NAME(zgeev)(&jobvl,&jobvr,&n,cpA.data(),&n,pd,pvl,&nevecl,pvr,&nevecr,work.data(),&lwork,rwork.data(),&info);
#endif
    return info;
    }

} //namespace itensor

