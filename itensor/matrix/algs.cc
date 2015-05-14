//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#include "algs.h"
#include "lapack_wrap.h"
#include <limits>
#include "count.h"
#include "slicemat.h"

using std::move;
using std::sqrt;

namespace itensor {

//
// diagSymmetric
//

void
diagSymmetric(const MatRefc& M,
              const MatRef& U,
              const VecRef& d)
    {
    LAPACK_INT N = M.Ncols();
    if(N < 1) throw std::runtime_error("diagSymmetric: 0 dimensional matrix");
    if(N != M.Nrows())
        {
        printfln("M is %dx%d",M.Nrows(),M.Ncols());
        throw std::runtime_error("diagSymmetric: Input Matrix must be square");
        }

#ifdef DEBUG
    if(!(U.Nrows()== N && U.Ncols() == N)) 
        throw std::runtime_error("diagSymmetric: U should have same dims as M");
    if(d.size() != N)
        throw std::runtime_error("diagSymmetric: d size should be linear size of M");
    if(!U.contiguous())
        throw std::runtime_error("diagSymmetric: U must be contiguous");
    if(!d.contiguous())
        throw std::runtime_error("diagSymmetric: d must be contiguous");
#endif

    //Set U = -M so eigenvalues will be sorted from largest to smallest
    auto pM = M.cbegin();
    for(auto& el : U) 
        { 
        el = -(*pM); 
        ++pM; 
        }

    LAPACK_INT info = 0;
    dsyev_wrapper('V','U',N,U.data(),d.data(),info);
    if(info != 0) throw std::runtime_error("Error condition in diagSymmetric");

    //Correct the signs of the eigenvalues:
    d *= -1;
    }

void
diagSymmetric(MatRefc M,
              Mat& U,
              Vec& d)
    {
    U.resize(M.Nrows(),M.Ncols());
    d.resize(M.Nrows());
    auto Uref = makeRef(U);
    auto dref = makeRef(d);
    diagSymmetric(M,Uref,dref);
    }

//
// orthog
//

void 
orthog(MatRef M, long num, long numpass)
    {
    if(num == -1) num = M.Ncols();
#ifdef DEBUG
    //if(num > M.Nrows() || (num == 0 && M.Ncols() > M.Nrows()))
    //    throw std::runtime_error("orthog: Ncols() > M.Nrows()");
#endif

    long nkeep = -1;// Orthogonalize to at most the column dim 
    if (num > 0 && num <= M.Ncols() && num <= M.Nrows())
        {
        nkeep = num;
        }
    else
        {
        nkeep = std::min(M.Nrows(), M.Ncols());
        }

    Vec dots(nkeep);
    MatRef Mcols;
    VecRef dotsref, 
           coli;
    for(auto i : count1(nkeep))
        {
        coli = column(M,i);
        auto nrm = norm(coli);
        if(nrm == 0.0)
            {
            randomize(coli);
            nrm = norm(coli);
            }
        coli /= nrm;
        if(i == 1) continue;

        Mcols = columns(M,1,i-1);
        dotsref = subVector(dots,1,i-1);
        for(auto pass : count1(numpass))
            {
            dotsref &= transpose(Mcols) * coli;
            coli -= Mcols * dotsref;
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

//
// SVD
//

//#define CHKSVD

void 
checksvd(const MatRefc& A, const MatRefc& U, const VecRefc& D, const MatRefc& V)
    {
    Mat Ach(U);
    for(long i = 1; i <= D.size(); ++i) column(Ach,i) *= D(i);
    Ach = Ach * transpose(V);
    Ach -= A;
    printfln("relative error with sqrt in low level svd is %.5E",norm(Ach)/norm(A));
    }

void
SVDRef(const MatRefc& M,
       const MatRef& U, 
       const VecRef& D, 
       const MatRef& V,
       Real thresh,
       Real northpass)
    {
    auto Mr = M.Nrows(), 
         Mc = M.Ncols();

    if(Mr > Mc)
        {
        SVDRef(transpose(M),V,D,U,thresh,northpass);
#ifdef CHKSVD
        checksvd(M,U,D,V);
#endif
        return;
        }

#ifdef DEBUG
    if(!(U.Nrows()==Mr && U.Ncols()==Mr)) 
        throw std::runtime_error("SVD (ref version), wrong size of U");
    if(!(V.Nrows()==Mc && V.Ncols()==Mr)) 
        throw std::runtime_error("SVD (ref version), wrong size of V");
    if(D.size()!=Mr)
        throw std::runtime_error("SVD (ref version), wrong size of D");
#endif

    //Form 'density matrix' rho
    Mat rho = M * transpose(M);

    //Diagonalize rho: evals are squares of singular vals
    diagSymmetric(rho,U,D);

    for(auto& el : D) // el = sqrt(fabs(el));
        {
        //This formula zeroes out any negative evals,
        //smaller error than taking their abs value
        el = sqrt((fabs(el)+el)/2);
        }

    //Put result of Mt*U==(V*D) in V storage
    mult(transpose(M),U,V);
    //Orthogonalize cols of V*D to obtain V
    orthog(V,Mr,northpass); //2 is the number of orthog passes

    long start = 2;
    auto D1t = D(1)*thresh;
    for(; start < Mr; ++start)
        {
        if(D(start) < D1t) break;
        }

    if(start >= (Mr-1)) 
        {
#ifdef CHKSVD
        checksvd(M,U,D,V);
#endif
        return;
        }

    //
    //Recursively SVD part of B 
    //for greater final accuracy
    //

    auto n = Mr-start+1;

    //reuse storage of rho to hold mv=M*columns(V,start,Mr)
    Mat mv = move(rho);
    mv.resize(Mr,n);
    mult(M,columns(V,start,Mr),mv);

    //b should be close to diagonal
    //but may not be perfect - fix it up below
    Mat b = rows(transpose(U),start,Mr)*mv;
   
    auto d = subVector(D,start,Mr);
    Mat u(n,n),
        v(n,n);
    SVDRef(b,u,d,v,thresh,northpass);

    auto Uu = move(mv);
    mult(columns(U,start,Mr),u,Uu);
    columns(U,start,Mr) &= Uu;

    columns(V,start,Mr) &= columns(V,start,Mr) * v;

#ifdef CHKSVD
	checksvd(M,U,D,V);
#endif

    return;
    }

void
SVD(const MatRefc& M,
    Mat& U, 
    Vec& D, 
    Mat& V,
    Real thresh,
    Real northpass)
    {
    auto Mr = M.Nrows(),
         Mc = M.Ncols();
    auto nsv = std::min(Mr,Mc);
    U.resize(Mr,nsv);
    V.resize(Mc,nsv);
    D.resize(nsv);
    SVDRef(M,U,D,V,thresh,northpass);
    }

} //namespace itensor
