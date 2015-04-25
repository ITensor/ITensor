//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "matrix_algs.h"
#include "lapack_wrap.h"
#include <limits>
#include "count.h"

namespace itensor {

void
diagSymmetric(const matrixref& M,
              matrixref& U,
              vecref& d)

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
    dsyev_wrapper('V','U',N,U.store(),d.store(),info);
    if(info != 0) throw std::runtime_error("Error condition in diagSymmetric");

    //Correct the signs of the eigenvalues:
    d *= -1;
    }

void
diagSymmetric(const matrixref& M,
              matrix& U,
              vec& d)
    {
    if(U.Nrows() != M.Nrows() || U.Ncols() != M.Ncols())
        U = matrix(M.Nrows(),M.Ncols());
    if(d.size() != M.Nrows()) d = vec(M.Nrows());
    matrixref& Uref = U;
    vecref& dref = d;
    diagSymmetric(M,Uref,dref);
    }

void 
orthog(const matrixref& M, long num, long numpass)
    {
    if(num == -1) num = M.Ncols();
#ifdef DEBUG
    if(M.readOnly()) throw std::runtime_error("orthog: matrixref is read only");
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

    vec dots(nkeep);
    matrixref Mcols;
    vecref dotsref, 
           coli;
    for(auto i : count1(nkeep))
        {
        coli = column(M,i);
        auto nrm = norm(coli);
        if(nrm == 0.0)
            {
            coli.randomize();
            nrm = norm(coli);
            }
        coli /= nrm;
        if(i == 1) continue;

        Mcols = columns(M,1,i-1);
        dotsref = subVector(dots,1,i-1);
        for(auto pass : count1(numpass))
            {
            dotsref = Mcols.t() * coli;
            coli -= Mcols * dotsref;
            auto nrm = norm(coli);
            if(nrm < 1E-3) --pass; //orthog is suspect
            if(nrm < 1E-10) // What if a subspace was zero in all vectors?
                {
                coli.randomize();
                nrm = norm(coli);
                }
            coli /= nrm;
            }
        }
    }

}; //namespace itensor
