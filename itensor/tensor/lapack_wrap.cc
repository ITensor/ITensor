#include "itensor/tensor/lapack_wrap.h"
#include "itensor/tensor/permutecplx.h"

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
#ifdef PLATFORM_macos
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
#ifdef PLATFORM_macos
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
void 
gemm_wrapper(bool transa, 
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
//#define USE_ZGEMM
#define METHOD1
//#define METHOD2

#ifdef USE_ZGEMM
    TIMER_START(31)
    auto palpha = (void*)(&alpha); 
    auto pbeta = (void*)(&beta); 
    cblas_zgemm(CblasColMajor,at,bt,m,n,k,palpha,(void*)A,lda,(void*)B,ldb,pbeta,(void*)C,m);
    TIMER_STOP(31)
#endif

#ifdef METHOD2
    if(not (alpha == Cplx(1.,0.) && beta == Cplx(0.,0.)))
        {
        throw std::runtime_error("alpha,beta must be 1,0");
        }

    auto aCsize = m*k;
    auto bCsize = k*n;
    auto cCsize = m*n;
    //WARNING: logically const, not thread safe!
    auto Ad = reinterpret_cast<Real*>(const_cast<Cplx*>(A));
    auto Bd = reinterpret_cast<Real*>(const_cast<Cplx*>(B));
    auto Cd = reinterpret_cast<Real*>(C);
    auto arb = Ad;
    auto aib = Ad+aCsize;
    auto brb = Bd;
    auto bib = Bd+bCsize;
    auto crb = Cd;
    auto cib = Cd+cCsize;
        //print("re(A):");
        //for(auto i = arb; i != aib; ++i)
        //    {
        //    print(" ",*i);
        //    }
        //println();
        //print("im(A):");
        //for(auto i = aib; i != aib+aCsize; ++i)
        //    {
        //    print(" ",*i);
        //    }
        //println();
    toCplx(Ad,aCsize,false);
    toCplx(Bd,bCsize,false);

    //print("re(A):");
    //for(auto i = arb; i != aib; ++i)
    //    {
    //    print(" ",*i);
    //    }
    //println();
    //print("im(A):");
    //for(auto i = aib; i != aib+aCsize; ++i)
    //    {
    //    print(" ",*i);
    //    }
    //println();

    cblas_dgemm(CblasColMajor,at,bt,m,n,k,+1,arb,lda,brb,ldb,0.,crb,m);
    cblas_dgemm(CblasColMajor,at,bt,m,n,k,-1,aib,lda,bib,ldb,1.,crb,m);
    cblas_dgemm(CblasColMajor,at,bt,m,n,k,+1,arb,lda,bib,ldb,0.,cib,m);
    cblas_dgemm(CblasColMajor,at,bt,m,n,k,+1,aib,lda,brb,ldb,1.,cib,m);
    toCplx(Ad,aCsize,true);
    toCplx(Bd,bCsize,true);
    toCplx(Cd,cCsize,true);
#endif //METHOD2

#ifdef METHOD1
    if(not (alpha == Cplx(1.,0.) && beta == Cplx(0.,0.)))
        {
        throw std::runtime_error("alpha,beta must be 1,0");
        }

    auto arsize = m*k;
    auto brsize = k*n;
    auto crsize = m*n;
    //TIMER_START(32)
    //auto a = std::vector<Real>(2*arsize);
    //auto b = std::vector<Real>(2*brsize);
    //auto c = std::vector<Real>(2*crsize);
    //Real* const arb = a.data();
    //Real* const aib = arb+arsize;
    //Real* const brb = b.data();
    //Real* const bib = brb+brsize;
    //Real* const crb = c.data();
    //Real* const cib = crb+crsize;

    //auto d = std::vector<Real>(2*(arsize+brsize+crsize));
    //Real* const arb = d.data();
    //Real* const aib = arb+arsize;
    //Real* const brb = aib+arsize;
    //Real* const bib = brb+brsize;
    //Real* const crb = bib+brsize;
    //Real* const cib = crb+crsize;

    auto d = std::vector<Real>(arsize+brsize+crsize);
    Real* const ab = d.data();
    Real* const ae = d.data()+arsize;
    Real* const bb = ae;
    Real* const be = bb+brsize;
    Real* const cb = be;
    Real* const ce = cb+crsize;
    //TIMER_STOP(32)

    struct Task
        {
        bool copyC = false;
        size_t Apart = 0,
               Bpart = 0,
               Cpart = 0;
        Real alpha = 0,
             beta  = 0;
        Task(size_t Ap,
             size_t Bp,
             size_t Cp,
             Real a,
             Real b)
          : Apart(Ap),
            Bpart(Bp),
            Cpart(Cp),
            alpha(a),
            beta(b)
            { }
        Task(size_t Cp)
           : copyC(true),Cpart(Cp)
            { }
        };
    //beta is zero because we're writing c to a buffer
    std::array<Task,6> tasks = {{Task(1,1,0,-1.,0.),
                                 Task(0,0,0,+1.,1.),
                                 Task(0),
                                 Task(1,0,1,+1.,0.),
                                 Task(0,1,1,+1.,1.),
                                 Task(1)
                                }};

    auto Ad = reinterpret_cast<const Real*>(A);
    auto Bd = reinterpret_cast<const Real*>(B);
    auto Cd = reinterpret_cast<Real*>(C);
    for(auto t : tasks)
        {
        if(t.copyC)
            {
            auto Cb = Cd+t.Cpart;
            for(auto cb_=cb; cb_ != ce; ++cb_,Cb+=2)
                {
                *Cb = *cb_;
                }
            }
        else
            {
            auto Ab = Ad+t.Apart;
            for(auto ab_=ab; ab_ != ae; ++ab_,Ab+=2)
                {
                *ab_ = *Ab;
                }
            auto Bb = Bd+t.Bpart;
            for(auto bb_=bb; bb_ != be; ++bb_,Bb+=2)
                {
                *bb_ = *Bb;
                }
            cblas_dgemm(CblasColMajor,at,bt,m,n,k,t.alpha,ab,lda,bb,ldb,t.beta,cb,m);
            }
        }

    //TIMER_START(35)
    //cblas_dgemm(CblasColMajor,at,bt,m,n,k,+1,arb,lda,brb,ldb,1.,crb,m);
    //cblas_dgemm(CblasColMajor,at,bt,m,n,k,+1,aib,lda,brb,ldb,0.,cib,m);
    //cblas_dgemm(CblasColMajor,at,bt,m,n,k,+1,arb,lda,bib,ldb,1.,cib,m);
    //TIMER_STOP(35)

    //TIMER_START(34)
    //auto Ab = A;
    //for(auto arb_=arb,aib_=aib; arb_ != aib; ++arb_,++aib_,++Ab)
    //    {
    //    *arb_ = reinterpret_cast<const Real*>(Ab)[0];
    //    *aib_ = reinterpret_cast<const Real*>(Ab)[1];
    //    }
    //auto Bb = B;
    //for(auto brb_=brb,bib_=bib; brb_ != bib; ++brb_,++bib_,++Bb)
    //    {
    //    *brb_ = reinterpret_cast<const Real*>(Bb)[0];
    //    *bib_  = reinterpret_cast<const Real*>(Bb)[1];
    //    }
    //TIMER_STOP(34)
    //cblas_dgemm(CblasColMajor,at,bt,m,n,k,-1,aib,lda,bib,ldb,0.,crb,m);
    //TIMER_START(34)
    //auto Cb = C;
    //for(auto crb_=crb,cib_=cib; crb_ != cib; ++crb_,++cib_,++Cb)
    //    {
    //    reinterpret_cast<Real*>(Cb)[0] = *crb_;
    //    reinterpret_cast<Real*>(Cb)[1] = *cib_;
    //    }
    //TIMER_STOP(34)

#endif//METHOD1

#else //use Fortran zgemm
    auto *npA = const_cast<Cplx*>(A);
    auto *npB = const_cast<Cplx*>(B);
    auto *pA = reinterpret_cast<LAPACK_COMPLEX*>(npA);
    auto *pB = reinterpret_cast<LAPACK_COMPLEX*>(npB);
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
    }

//
// dgemv - matrix*vector multiply
//
void 
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
void 
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
void 
dscal_wrapper(LAPACK_INT N,
              LAPACK_REAL alpha,
              LAPACK_REAL* data,
              LAPACK_INT inc)
    {
#ifdef PLATFORM_macos
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
LAPACK_INT 
zheev_wrapper(LAPACK_INT N,        //number of cols of A
              LAPACK_COMPLEX *A,    //matrix A, on return contains eigenvectors
              LAPACK_REAL *d)       //eigenvalues on return
    {
    char jobz = 'V';
    char uplo = 'U';
    LAPACK_INT lwork = std::max(1,3*N-1);//max(1, 1+6*N+2*N*N);
    std::vector<LAPACK_COMPLEX> work(lwork);
    std::vector<LAPACK_REAL> rwork(lwork);
    LAPACK_INT info = 0;
#ifdef PLATFORM_acml
    LAPACK_INT jobz_len = 1;
    LAPACK_INT uplo_len = 1;
    F77NAME(zheev)(&jobz,&uplo,&N,A,&N,d,work.data(),&lwork,rwork.data(),&info,jobz_len,uplo_len);
#else
    F77NAME(zheev)(&jobz,&uplo,&N,A,&N,d,work.data(),&lwork,rwork.data(),&info);
#endif
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
void 
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
void 
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

} //namespace itensor

