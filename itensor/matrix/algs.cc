//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#include "algs.h"
#include "lapack_wrap.h"
#include <limits>
#include "count.h"
#include "slicemat.h"

namespace itensor {

//
// diagSymmetric
//

void
diagSymmetric(MatRefc M,
              MatRef U,
              VecRef d)
    {
    LAPACK_INT N = M.Ncols();
    if(N < 1) throw std::runtime_error("diagSymmetric: 0 dimensional matrix");
    if(N != M.Nrows())
        {
        printfln("M is %dx%d",M.Nrows(),M.Ncols());
        throw std::runtime_error("diagSymmetric: Input Matrix must be square");
        }

#ifdef DEBUG
    if(U.Nrows() != N || U.Ncols() != N) 
        throw std::runtime_error("diagSymmetric: U should have same dims as M");
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

    LAPACK_INT info;
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
    if(!(U.Nrows()==M.Nrows() && U.Ncols()==M.Ncols()))
        {
        U = Mat(M.Nrows(),M.Ncols());
        }
    if(d.size() != M.Nrows()) 
        {
        d = Vec(M.Nrows());
        }
    diagSymmetric(M,makeRef(U),makeRef(d));
    }

//
// orthog
//

void 
orthog(MatRef M, long num, long numpass)
    {
    if(num == -1) num = M.Ncols();
#ifdef DEBUG
    if(num > M.Nrows() || (num == 0 && M.Ncols() > M.Nrows()))
        throw std::runtime_error("orthog: Ncols() > M.Nrows()");
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
            auto nrm = norm(coli);
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

#define CHKSVD

/*
void 
checksvd(const MatRefc& A, const MatRefc& U, const VecRefc& D, const MatRefc& V)
    {
    Mat Ach(U);
    for(long i = 1; i <= D.size(); ++i) column(Ach,i) *= D(i);
    Ach = Ach * V;
    println("A = \n",A);
    println("Ach = \n",Ach);
    Ach -= A;
    printfln("relative error with sqrt in low level svd is %.5E",norm(Ach)/norm(A));
    }

void
SVD(MatRefc A,
    MatRef U, 
    VecRef D, 
    MatRef V,
    Real thresh)
    {
    println("In SVD");
    auto n = A.Nrows(), 
         m = A.Ncols();

    if(n > m)
        {
        auto At = transpose(A);
        auto Ut = transpose(U);
        auto Vt = transpose(V);
        SVD(At,Vt,D,Ut,thresh);
#ifdef CHKSVD
        checksvd(A,U,D,V);
#endif
        return;
        }

    //Form 'density matrix' rho
    Mat rho = A * transpose(A);

    Vec evals(n);
    diagSymmetric(rho,U,evals);

    //Form Vt and fix up its orthogonality
    //(Vt is transpose of V)

    V &= transpose(U) * A;
    orthog(V,n,2); //2 is the number of orthog passes
    //V &= transpose(Vt);

    Mat Vt = transpose(A) * U;
    orthog(Vt,n,2); //2 is the number of orthog passes
    V &= transpose(Vt);

    //B should be close to diagonal
    //but may not be perfect - fix
    //it up below
    const Mat B = transpose(U) * A * Vt;

    D &= diagonal(B);

    if(D(1) == 0 || thresh == 0)
        {
#ifdef CHKSVD
        checksvd(A,U,D,V);
#endif
        return;
        }

    long start = 2;
    auto D1 = D(1);
    for(; start < n; ++start)
        {
        if(D(start)/D1 < thresh) break;
        }

    if(start >= (n-1)) 
        {
#ifdef CHKSVD
        checksvd(A,U,D,V);
#endif
        return;
        }

    //
    //Recursively SVD part of B 
    //for greater final accuracy
    //

    auto b = subMatrix(B,start,n,start,n);

    Mat u,
        v;
    Vec d;
    SVD(b,u,d,v,thresh);

    subVector(D,start,n) &= d;

    subMatrix(U,1,n,start,n) &= subMatrix(U,1,n,start,n) * u;

    subMatrix(V,start,n,1,m) &= v * subMatrix(V,start,n,1,m);

#ifdef CHKSVD
	checksvd(A,U,D,V);
#endif

    return;
    }

void
SVD(MatRefc A,
    Mat& U, 
    Vec& D, 
    Mat& V,
    Real thresh)
    {
    auto nsv = std::min(A.Nrows(),A.Ncols());
    if(!(U.Nrows()==A.Nrows() && U.Ncols()==nsv)) 
        {
        U = Mat(A.Nrows(),nsv);
        }
    if(!(V.Nrows()==nsv && V.Ncols()==A.Ncols())) 
        {
        V = Mat(nsv,A.Ncols());
        }
    if(D.size() != nsv) 
        {
        D = Vec(nsv);
        }
    SVD(A,MatRef(U),VecRef(D),MatRef(V),thresh);
    }
*/

}; //namespace itensor
