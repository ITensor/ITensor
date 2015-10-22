//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#include <limits>
#include <stdexcept>
#include <tuple>
#include "itensor/tensor/lapack_wrap.h"
#include "itensor/tensor/algs.h"
#include "itensor/util/count.h"

using std::move;
using std::sqrt;
using std::tuple;
using std::make_tuple;
using std::tie;

namespace itensor {

//Helper for diagSymmetric
template<typename Iter>
void
copyNegElts(Iter m,
            MatrixRef const& U)
    {
    auto ue = U.data()+U.size();
    for(auto u = U.data(); u != ue; ++u,++m)
        {
        *u = -(*m);
        }
    }

//Helper for diagHermitian
template<typename Iter>
void
copyNegElts(Iter mre,
            Iter mim,
            std::vector<LAPACK_COMPLEX> & Mc)
    {
    for(auto& z : Mc)
        {
        realRef(z) = -(*mre);
        imagRef(z) = -(*mim);
        ++mre;
        ++mim;
        }
    }

//
// diagSymmetric
//

void
diagSymmetric(MatrixRefc const& M,
              MatrixRef  const& U,
              VectorRef  const& d)
    {
    auto N = ncols(M);
    if(N < 1) throw std::runtime_error("diagSymmetric: 0 dimensional matrix");
    if(N != nrows(M))
        {
        printfln("M is %dx%d",nrows(M),ncols(M));
        throw std::runtime_error("diagSymmetric: Input Matrix must be square");
        }

#ifdef DEBUG
    if(!(nrows(U) == N && ncols(U) == N)) 
        throw std::runtime_error("diagSymmetric: U should have same dims as M");
    if(d.size() != N)
        throw std::runtime_error("diagSymmetric: d size should be linear size of M");
    if(!isContiguous(U))
        throw std::runtime_error("diagSymmetric: U must be contiguous");
    if(!isContiguous(d))
        throw std::runtime_error("diagSymmetric: d must be contiguous");
#endif

    //Set U = -M so eigenvalues will be sorted from largest to smallest
    if(isContiguous(M)) copyNegElts(M.data(),U);
    else                copyNegElts(M.cbegin(),U);

    LAPACK_INT info = 0;
    dsyev_wrapper('V','U',N,U.data(),d.data(),info);
    if(info != 0) 
        {
        println("M = \n",M);
        throw std::runtime_error("Error condition in diagSymmetric");
        }

    //Correct the signs of the eigenvalues:
    d *= -1;
    }

void
diagSymmetric(MatrixRefc const& M,
              Matrix          & U,
              VectorRef  const& d)
    {
    resize(U,nrows(M),ncols(M));
    diagSymmetric(M,makeRef(U),d);
    }

void
diagSymmetric(MatrixRefc const& M,
              Matrix          & U,
              Vector          & d)
    {
    resize(U,nrows(M),ncols(M));
    resize(d,nrows(M));
    diagSymmetric(M,makeRef(U),makeRef(d));
    }

void
diagHermitian(MatrixRefc const& Mre,
              MatrixRefc const& Mim,
              MatrixRef  const& Ure,
              MatrixRef  const& Uim,
              VectorRef  const& d)
    {
    auto N = ncols(Mre);
    if(N != nrows(Mre))
        {
        printfln("Mre is %dx%d",nrows(Mre),ncols(Mre));
        throw std::runtime_error("diagHermitian: Input Matrix must be square");
        }
    if(N != nrows(Mim) || N != ncols(Mim))
        {
        printfln("Mim is %dx%d",nrows(Mim),ncols(Mim));
        throw std::runtime_error("diagHermitian: Input Matrix must be square, and real and imag part same size");
        }

#ifdef DEBUG
    if(N < 1) throw std::runtime_error("diagHermitian: 0 dimensional matrix");
    if(!(nrows(Ure) == N && ncols(Ure) == N)) 
        throw std::runtime_error("diagHermitian: Ure should have same dims as M");
    if(!(nrows(Uim) == N && ncols(Uim) == N)) 
        throw std::runtime_error("diagHermitian: Uim should have same dims as M");
    if(d.size() != N)
        throw std::runtime_error("diagHermitian: d size should be linear size of M");
    if(!isContiguous(Ure))
        throw std::runtime_error("diagHermitian: Ure must be contiguous");
    if(!isContiguous(Uim))
        throw std::runtime_error("diagHermitian: Uim must be contiguous");
    if(!isContiguous(d))
        throw std::runtime_error("diagHermitian: d must be contiguous");
#endif


    //Set Mc = -M so eigenvalues will be sorted from largest to smallest
    auto Mc = std::vector<LAPACK_COMPLEX>(N*N);
    if(isContiguous(Mre) && isContiguous(Mim))
        {
        copyNegElts(Mre.data(),Mim.data(),Mc);
        }
    else
        {
        copyNegElts(Mre.cbegin(),Mim.cbegin(),Mc);
        }

    auto info = zheev_wrapper(N,Mc.data(),d.data());
    if(info != 0) 
        {
        throw std::runtime_error("Error condition in diagHermitian");
        }

    //Correct eigenvalue signs
    d *= -1;

    //Following code assumes Ure and Uim are contiguous
    auto ur = Ure.data();
    auto ui = Uim.data();
    for(auto& z : Mc)
        {
        (*ur) = realRef(z);
        (*ui) = imagRef(z);
        ++ur;
        ++ui;
        }
    }

void
diagHermitian(MatrixRefc const& Mre,
              MatrixRefc const& Mim,
              Matrix          & Ure,
              Matrix          & Uim,
              VectorRef  const& d)
    {
    resize(Ure,nrows(Mre),ncols(Mre));
    resize(Uim,nrows(Mre),ncols(Mre));
    diagHermitian(Mre,Mim,makeRef(Ure),makeRef(Uim),d);
    }

void
diagHermitian(MatrixRefc const& Mre,
              MatrixRefc const& Mim,
              Matrix          & Ure,
              Matrix          & Uim,
              Vector          & d)
    {
    resize(Ure,nrows(Mre),ncols(Mre));
    resize(Uim,nrows(Mre),ncols(Mre));
    resize(d,nrows(Mre));
    diagHermitian(Mre,Mim,makeRef(Ure),makeRef(Uim),makeRef(d));
    }

//
// orthog
//

void 
orthog(MatrixRef M, 
       size_t numpass)
    {
    auto nkeep = std::min(nrows(M), ncols(M));
    auto dots = Vector(nkeep);
    for(auto i : count(nkeep))
        {
        //normalize column i
        auto coli = column(M,i);
        auto nrm = norm(coli);
        if(nrm == 0.0)
            {
            randomize(coli);
            nrm = norm(coli);
            }
        coli /= nrm;
        if(i == 0) continue;

        auto Mcols = columns(M,0,i);
        auto dotsref = subVector(dots,0,i);
        for(auto pass : count1(numpass))
            {
            // does dotsref &= transpose(Mcols) * coli:
            mult(transpose(Mcols),coli,dotsref);
            // does coli -= Mcols * dotsref:
            multSub(Mcols,dotsref,coli);
            nrm = norm(coli);
            if(nrm < 1E-3) --pass; //orthog is suspect
            if(nrm < 1E-10) // What if a subspace was zero in all vectors?
                {
                randomize(coli);
                nrm = norm(coli);
                }
            coli /= nrm;
            }
        }
    }

Real static
sqr(Real x) { return x*x; }

void 
orthog(MatrixRef Mr,
       MatrixRef Mi,
       size_t numpass)
    {
    auto nkeep = std::min(nrows(Mr), ncols(Mr));
    auto Dr = Vector(nkeep);
    auto Di = Vector(nkeep);
    auto cnorm = [](VectorRefc const& r,
                    VectorRefc const& i)
        {
        return sqrt(sqr(norm(r))+sqr(norm(i)));
        };
    for(auto n : count(nkeep))
        {
        //normalize column n
        auto cr = column(Mr,n);
        auto ci = column(Mi,n);
        auto nrm = cnorm(cr,ci);
        if(nrm == 0.0)
            {
            randomize(cr);
            randomize(ci);
            nrm = cnorm(cr,ci);
            }
        cr /= nrm;
        ci /= nrm;
        if(n == 0) continue;

        auto mr = columns(Mr,0,n);
        auto mi = columns(Mi,0,n);
        auto dr = subVector(Dr,0,n);
        auto di = subVector(Di,0,n);
        for(auto pass : count1(numpass))
            {
            //// does dotsref &= transpose(Mcols) * coli:
            //mult(transpose(Mcols),coli,dotsref);
            dr &= transpose(mr)*cr+transpose(mi)*ci;
            di &= transpose(mr)*ci-transpose(mi)*cr;
            cr -= mr*dr-mi*di;
            ci -= mr*di+mi*dr;

            nrm = cnorm(cr,ci);
            if(nrm < 1E-3) --pass; //orthog is suspect
            if(nrm < 1E-10) // What if a subspace was zero in all vectors?
                {
                randomize(cr);
                randomize(ci);
                nrm = cnorm(cr,ci);
                }
            cr /= nrm;
            ci /= nrm;
            }
        }
    }

//
// SVD
//

//#define CHKSVD

void 
checksvd(MatrixRefc const& A, 
         MatrixRefc const& U, 
         VectorRefc const& D, 
         MatrixRefc const& V)
    {
    Matrix Ach(U);
    for(auto i : count1(D.size())) column(Ach,i) *= D(i);
    Ach = Ach * transpose(V);
    Ach -= A;
    printfln("relative error with sqrt in low level svd is %.5E",norm(Ach)/norm(A));
    }

tuple<bool,size_t> // == (done, start)
checkSVDDone(VectorRefc const& D,
             Real thresh)
    {
    auto N = D.size();
    if(N <= 1 || thresh <= 0) 
        {
        println("Got zero thresh");
        return make_tuple(true,1);
        }
    auto D1t = D(0)*thresh;
    size_t start = 1;
    for(; start < N; ++start)
        {
        if(D(start) < D1t) break;
        }

    if(start >= (N-1)) 
        return make_tuple(true,start);

    return make_tuple(false,start);
    }

void
SVDRefImpl(MatrixRefc const& M,
           MatrixRef  const& U, 
           VectorRef  const& D, 
           MatrixRef  const& V,
           Real thresh,
           int depth = 0)
    {
    auto Mr = nrows(M), 
         Mc = ncols(M);

    if(Mr > Mc)
        {
        SVDRefImpl(transpose(M),V,D,U,thresh,depth);
#ifdef CHKSVD
        checksvd(M,U,D,V);
#endif
        return;
        }

#ifdef DEBUG
    if(!(nrows(U)==Mr && ncols(U)==Mr)) 
        throw std::runtime_error("SVD (ref version), wrong size of U");
    if(!(nrows(V)==Mc && ncols(V)==Mr)) 
        throw std::runtime_error("SVD (ref version), wrong size of V");
    if(D.size()!=Mr)
        throw std::runtime_error("SVD (ref version), wrong size of D");
#endif

    //Form 'density matrix' rho
    auto rho = M * transpose(M);

    //Diagonalize rho: evals are squares of singular vals
    diagSymmetric(rho,U,D);

    for(auto& el : D)
        {
        if(el < 0) el = 0.;
        else       el = std::sqrt(el);
        }

    size_t nlarge = 0;
    auto rthresh = D(0)*thresh;
    for(decltype(Mr) n = 0; n < Mr; ++n)
        {
        if(D(n) < rthresh)
            {
            nlarge = n;
            break;
            }
        }

    //Put result of Mt*U==(V*D) in V storage
    mult(transpose(M),U,V);

    for(decltype(nlarge) n = 0; n < nlarge; ++n)
        {
        column(V,n) /= D(n);
        }

    if(nlarge < Mr)
        {
        //Much more accurate than dividing
        //by smallest singular values
        orthog(columns(V,nlarge,Mr),2);
        }

    bool done = false;
    size_t start = 1;
    tie(done,start) = checkSVDDone(D,thresh);

    if(done) return;

    //
    //Recursively SVD part of B 
    //for greater final accuracy
    //

    auto n = Mr-start;

     //   {
     //   //println("Method 1");
     //   //TEST VERSION - SLOW!
     //   auto B = transpose(U)*M*V;
     //   auto b = Matrix{subMatrix(B,start,Mr,start,Mr)};

     //   auto d = subVector(D,start,Mr);
     //   Matrix u(n,n),
     //          v(n,n);
     //   SVDRefImpl(b,u,d,v,thresh,1+depth);

     //   auto ns = d.size();
     //   auto dd = Matrix(ns,ns);
     //   diagonal(dd) &= d;

     //   auto nu = columns(U,start,Mr);
     //   auto tmpu = nu * u;
     //   nu &= tmpu;

     //   auto nv = columns(V,start,Mr);
     //   auto tmpv = nv * v;
     //   nv &= tmpv;
     //   }

    //reuse storage of rho to hold mv=M*columns(V,start,Mr)
    auto mv = move(rho);
    reduceCols(mv,n);

    auto u = columns(U,start,Mr);
    auto v = columns(V,start,Mr);

    //b should be close to diagonal
    //but may not be perfect - fix it up below
    mult(M,v,mv);
    auto b = transpose(u)*mv;
   
    auto d = subVector(D,start,Mr);
    Matrix bu(n,n),
           bv(n,n);
    SVDRef(b,bu,d,bv,thresh);

    //reuse storage to avoid allocations
    auto W = move(mv);
    mult(u,bu,W);
    u &= W;

    mult(v,bv,W);
    v &= W;

#ifdef CHKSVD
	checksvd(M,U,D,V);
#endif

    return;
    }

void
SVDRef(MatrixRefc const& M,
       MatrixRef  const& U, 
       VectorRef  const& D, 
       MatrixRef  const& V,
       Real thresh)
    {
    SVDRefImpl(M,U,D,V,thresh);
    }


void
SVD(MatrixRefc const& M,
    Matrix & U, 
    Vector & D, 
    Matrix & V,
    Real thresh)
    {
    auto Mr = nrows(M),
         Mc = ncols(M);
    auto nsv = std::min(Mr,Mc);
    resize(U,Mr,nsv);
    resize(V,Mc,nsv);
    resize(D,nsv);
    SVDRef(M,U,D,V,thresh);
    }

void
SVDRef(MatrixRefc const& Mre,
       MatrixRefc const& Mim,
       MatrixRef  const& Ure, 
       MatrixRef  const& Uim, 
       VectorRef  const& D, 
       MatrixRef  const& Vre,
       MatrixRef  const& Vim,
       Real thresh)
    {
    auto Mr = nrows(Mre), 
         Mc = ncols(Mim);

    if(Mr > Mc)
        {
        SVDRef(transpose(Mre),transpose(Mim),Vre,Vim,D,Ure,Uim,thresh);
        Uim *= -1;
        Vim *= -1;
        return;
        }

#ifdef DEBUG
    if(!(nrows(Mim)==Mr && ncols(Mim)==Mc)) 
        throw std::runtime_error("SVD (ref version), Mim must have same dims as Mre");
    if(!(nrows(Ure)==Mr && ncols(Ure)==Mr)) 
        throw std::runtime_error("SVD (ref version), wrong size of Ure");
    if(!(nrows(Uim)==Mr && ncols(Uim)==Mr)) 
        throw std::runtime_error("SVD (ref version), wrong size of Uim");
    if(!(nrows(Vre)==Mc && ncols(Vre)==Mr)) 
        throw std::runtime_error("SVD (ref version), wrong size of Vre");
    if(!(nrows(Vim)==Mc && ncols(Vim)==Mr)) 
        throw std::runtime_error("SVD (ref version), wrong size of Vim");
    if(D.size()!=Mr)
        throw std::runtime_error("SVD (ref version), wrong size of D");
#endif

    //Form 'density matrix' rho
    auto rhore = Mre*transpose(Mre) + Mim*transpose(Mim);
    auto rhoim = Mim*transpose(Mre) - Mre*transpose(Mim);

    //Diagonalize rho: evals are squares of singular vals
    diagHermitian(rhore,rhoim,Ure,Uim,D);

    for(auto& el : D)
        {
        if(el < 0) el = 0.;
        else       el = std::sqrt(el);
        }
    size_t nlarge = 0;
    auto rthresh = D(0)*thresh;
    for(decltype(Mr) n = 0; n < Mr; ++n)
        {
        if(D(n) < rthresh)
            {
            nlarge = n;
            break;
            }
        }

    //Compute Mt*U = V*D
    Vre &= transpose(Mre)*Ure + transpose(Mim)*Uim;
    Vim &= transpose(Mre)*Uim - transpose(Mim)*Ure;

    for(decltype(nlarge) n = 0; n < nlarge; ++n)
        {
        column(Vre,n) /= D(n);
        column(Vim,n) /= D(n);
        }
    if(nlarge < Mr)
        {
        //Much more accurate than dividing
        //by smallest singular values
        auto Vcr = columns(Vre,nlarge,Mr);
        auto Vci = columns(Vim,nlarge,Mr);
        orthog(Vcr,Vci,2);
        }

    bool done = false;
    size_t start = 1;
    tie(done,start) = checkSVDDone(D,thresh);
    if(done) return;

    //
    //Recursively SVD part of B 
    //for greater final accuracy
    //
    auto n = Mr-start;

    //{
    //println("Method 1");
    ////TEST VERSION - SLOW!
    //auto Tre = Mre*Vre - Mim*Vim;
    //auto Tim = Mre*Vim + Mim*Vre;
    //auto Bre = transpose(Ure)*Tre + transpose(Uim)*Tim;
    //auto Bim = transpose(Ure)*Tim - transpose(Uim)*Tre;

    //auto bre = Matrix{subMatrix(Bre,start,Mr,start,Mr)};
    //auto bim = Matrix{subMatrix(Bim,start,Mr,start,Mr)};

    //auto d = subVector(D,start,Mr);
    //Matrix ure(n,n),
    //       uim(n,n),
    //       vre(n,n),
    //       vim(n,n);
    //SVDRef(bre,bim,ure,uim,d,vre,vim,thresh);

    //auto nure = columns(Ure,start,Mr);
    //auto nuim = columns(Uim,start,Mr);
    //auto tmpre = nure*ure - nuim*uim;
    //auto tmpim = nuim*ure + nure*uim;
    //nure &= tmpre;
    //nuim &= tmpim;

    //auto nvre = columns(Vre,start,Mr);
    //auto nvim = columns(Vim,start,Mr);
    //tmpre = nvre*vre - nvim*vim;
    //tmpim = nvim*vre + nvre*vim;
    //nvre &= tmpre;
    //nvim &= tmpim;
    //}

    //reuse storage of rho to hold mv=M*columns(V,start,Mr)
    auto mvre = move(rhore);
    auto mvim = move(rhoim);
    reduceCols(mvre,n);
    reduceCols(mvim,n);

    auto ure = columns(Ure,start,Mr);
    auto uim = columns(Uim,start,Mr);
    auto vre = columns(Vre,start,Mr);
    auto vim = columns(Vim,start,Mr);

    mvre = Mre*vre - Mim*vim;
    mvim = Mre*vim + Mim*vre;

    auto utre = transpose(ure);
    auto utim = transpose(uim);

    //b (=ut*M*v) should be close to diagonal
    //but may not be perfect - fix it up below
    auto bre = utre*mvre + utim*mvim;
    auto bim = utre*mvim - utim*mvre;
    auto d = subVector(D,start,Mr);
    Matrix bure(n,n),
           buim(n,n),
           bvre(n,n),
           bvim(n,n);
    SVDRef(bre,bim,bure,buim,d,bvre,bvim,thresh);

    auto Nure = ure*bure-uim*buim;
    auto Nuim = ure*buim+uim*bure;
    ure &= Nure;
    uim &= Nuim;

    auto Nvre = vre*bvre-vim*bvim;
    auto Nvim = vre*bvim+vim*bvre;
    vre &= Nvre;
    vim &= Nvim;

#ifdef CHKSVD
	checksvd(M,U,D,V);
#endif

    return;
    }

void
SVD(MatrixRefc const& Mre,
    MatrixRefc const& Mim,
    Matrix & Ure, 
    Matrix & Uim, 
    Vector & D, 
    Matrix & Vre,
    Matrix & Vim,
    Real thresh)
    {
    auto Mr = nrows(Mre),
         Mc = ncols(Mim);
    auto nsv = std::min(Mr,Mc);
    resize(Ure,Mr,nsv);
    resize(Uim,Mr,nsv);
    resize(Vre,Mc,nsv);
    resize(Vim,Mc,nsv);
    resize(D,nsv);
    SVDRef(Mre,Mim,Ure,Uim,D,Vre,Vim,thresh);
    }


} //namespace itensor
