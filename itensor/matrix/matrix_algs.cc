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

    char jobz = 'V';
    char uplo = 'U';
    LAPACK_INT info;

    std::copy(M.cbegin(),M.cend(),U.begin());
    
    dsyev_wrapper(&jobz,&uplo,&N,U.store(),&N,d.store(),&info);

    if(info != 0) throw std::runtime_error("Error condition in diagSymmetric");

    //Transpose U before return
    U.applyTrans();
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
//#ifdef DEBUG
//    if(M.readOnly()) throw std::runtime_error("orthog: matrixref is read only");
//    if(num > M.Nrows() || (num == 0 && M.Ncols() > M.Nrows()))
//        throw std::runtime_error("orthog: Ncols() > M.Nrows()");
//#endif
//
//    long nkeep = -1;// Orthogonalize to at most the column dim 
//    if (num > 0 && num <= M.Ncols() && num <= M.Nrows())
//        {
//        nkeep = num;
//        }
//    else
//        {
//        nkeep = std::min(M.Nrows(), M.Ncols());
//        }
//
//    vec dots(nkeep);
//    matrixref Mcols;
//    vecref dotsref, 
//           coli;
//    for(auto i : count1(nkeep))
//        {
//        coli = column(M,i);
//        auto norm = norm(coli);
//        if(norm == 0.0)
//            {
//            for(auto& el : coli) el = detail::quickran();
//            norm = norm(coli);
//            }
//        coli /= norm;
//        if(i == 1) continue;
//
//        Mcols = columns(M,1,i-1);
//        dotsref = subVector(dots,1,i-1);
//        for(auto pass : count1(numpass))
//            {
//            dotsref = Mcols.t() * coli;
//            coli -= Mcols * dotsref;
//            auto norm = norm(coli);
//            if(norm < 1E-3)   // orthog is suspect
//            --pass;
//            if(norm < 1E-10)  // What if a subspace was zero in all vectors?
//                {
//                for(auto& el : coli) el = detail::quickran();
//                norm = norm(coli);
//                }
//            coli /= norm;
//            }
//        }
    }

}; //namespace itensor
