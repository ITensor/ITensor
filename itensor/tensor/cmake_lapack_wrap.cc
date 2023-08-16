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
#include <blas.hh>   // BLASPP
#include <lapack.hh> // LAPACKPP

#include "itensor/tensor/lapack_wrap.h"

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
    blas::axpy(n, alpha, X, incx, Y, incy);
    }

//
// dnrm2
//
LAPACK_REAL 
dnrm2_wrapper(LAPACK_INT N,
              const LAPACK_REAL* X,
              LAPACK_INT incx)
    {
      return blas::nrm2(N, X, incx);
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
      return blas::dot(N, X, incx, Y, incy);
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
      return blas::dot(N, X, incx, Y, incy);
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
    auto at = blas::Op::NoTrans,
    bt = blas::Op::NoTrans;
    if(transa) {
        at = blas::Op::Trans;
        lda = k;
    } if(transb){
        bt = blas::Op::Trans;
        ldb = n;
    }
    blas::gemm(blas::Layout::ColMajor, at,bt,m,n,k,alpha,A,lda,B,ldb,beta,C,m);
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
        auto at = blas::Op::NoTrans,
                bt = blas::Op::NoTrans;
        if(transa) {
            at = blas::Op::Trans;
            lda = k;
        } if(transb){
            bt = blas::Op::Trans;
            ldb = n;
        }
        blas::gemm(blas::Layout::ColMajor, at,bt,m,n,k,alpha,A,lda,B,ldb,beta,C,m);
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
    auto Tr = trans ? blas::Op::Trans : blas::Op::NoTrans;
    blas::gemv(blas::Layout::ColMajor, Tr, m, n, alpha, A, m, x, incx, beta, y, incy);
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
        auto Tr = trans ? blas::Op::Trans : blas::Op::NoTrans;
        blas::gemv(blas::Layout::ColMajor, Tr, m, n, alpha, A, m, x, incx, beta, y, incy);
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
    LAPACK_INT lda = n;
    info = lapack::syev(lapack::char2job(jobz), blas::char2uplo(uplo), n, A, lda, eigs);
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
    blas::scal(N, alpha, data, inc);
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
    LAPACK_INT l = std::min(*m,*n);
    *info = lapack::gesdd(lapack::char2job(*jobz), *m, *n, A, *m, s, u, *m, vt, l);
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
        LAPACK_INT l = std::min(*m,*n);
        *info = lapack::gesdd(lapack::char2job(*jobz), *m, *n, A, *m, s, u, *m, vt, l);
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
    LAPACK_INT l = std::min(*m,*n);
    *info = lapack::gesvd(lapack::char2job(*jobz), lapack::char2job(*jobz), *m, *n, A, *m, s, u, *m, vt, l);
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
        LAPACK_INT l = std::min(*m,*n);
        *info = lapack::gesvd(lapack::char2job(*jobz), lapack::char2job(*jobz), *m, *n, A, *m, s, u, *m, vt, l);
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
        *info = lapack::geqrf(*m, *n, A, *lda, tau);
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
        *info = lapack::orgqr(*m, *n, *k, A, *lda, tau);
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
    static_assert(sizeof(LAPACK_COMPLEX)==sizeof(Cplx),"LAPACK_COMPLEX and itensor::Cplx have different size");
    auto ptau = reinterpret_cast<Cplx*>(tau);
    *info = lapack::geqrf(*m, *n, A, *lda, ptau);
    tau = reinterpret_cast<LAPACK_COMPLEX*>(ptau);
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
    static_assert(sizeof(LAPACK_COMPLEX)==sizeof(Cplx),"LAPACK_COMPLEX and itensor::Cplx have different size");
    auto ptau = reinterpret_cast<Cplx *>(tau);
    *info = lapack::ungqr(*m, *n, *k, A, *lda, ptau);
    tau = reinterpret_cast<LAPACK_COMPLEX*>(ptau);
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
	std::vector<int64_t> ipiv(n);
	LAPACK_INT ldb = n;
    return lapack::gesv(n, nrhs, a, lda, ipiv.data(), b, ldb);
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
	LAPACK_INT lda = n;
	std::vector<int64_t> ipiv(n);
	LAPACK_INT ldb = n;
	auto info = lapack::gesv(n, nrhs, a, lda, ipiv.data(), b, ldb);
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
	return lapack::lange(lapack::char2norm(norm), m, n, a, m);
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
    return lapack::lange(lapack::char2norm(norm), m, n, a, m);
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
    char jobz = 'V';
    char uplo = 'U';
    lapack_int info = lapack::heev(lapack::char2job(jobz), blas::char2uplo(uplo), N, A, N, d);

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
    LAPACK_INT itype = 1;
    *info = lapack::sygv(itype, lapack::char2job(*jobz), blas::char2uplo(*uplo), *n, A, *n, B, *n, d);
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
    std::vector<Cplx> W(n);
    info = lapack::geev(lapack::char2job(jobvl), lapack::char2job(jobvr), n, cpA.data(), n, W.data(), vl, nevecl, vr, nevecr);
    auto ptr = W.data();
    for(size_t i = 0; i < n; ++i){
        *(dr + i) = std::real(*(ptr + i));
        *(di + i) = std::imag(*(ptr + i));
    }
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
    std::vector<Cplx> cpA;
    LAPACK_INT nevecl = (jobvl == 'V' ? n : 1);
    LAPACK_INT nevecr = (jobvr == 'V' ? n : 1);

    //Copy A data into cpA
    cpA.resize(n*n);
    auto pA = reinterpret_cast<Cplx const*>(A);
    std::copy(pA,pA+n*n,cpA.data());

    LAPACK_INT info = 0;
    info = lapack::geev(lapack::char2job(jobvl), lapack::char2job(jobvr), n,
                        cpA.data(), n, d, vl, nevecl, vr, nevecr);

    return info;
    }

} //namespace itensor

