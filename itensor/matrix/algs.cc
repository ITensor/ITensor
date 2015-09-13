//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#include <limits>
#include <stdexcept>
#include "itensor/matrix/lapack_wrap.h"
#include "itensor/matrix/algs.h"
#include "itensor/util/count.h"

using std::move;
using std::sqrt;

namespace itensor {

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
    if(isContiguous(M) && isContiguous(U))
        {
        daxpy_wrapper(M.size(),-1,M.data(),1,U.data(),1);
        }
    else
        {
        auto pM = M.cbegin();
        for(auto& el : U) 
            { 
            el = -(*pM); 
            ++pM; 
            }
        }

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
diagSymmetric(MatrixRefc M,
              Matrix & U,
              Vector & d)
    {
    resize(U,nrows(M),ncols(M));
    resize(d,nrows(M));
    auto Uref = makeRef(U);
    auto dref = makeRef(d);
    diagSymmetric(M,Uref,dref);
    }

//
// orthog
//

void 
orthog(MatrixRef M, size_t num, size_t numpass)
    {
    if(num == 0) num = ncols(M);
#ifdef DEBUG
    //if(num > nrows(M) || (num == 0 && ncols(M) > nrows(M)))
    //    throw std::runtime_error("orthog: ncols() > nrows()");
#endif

    size_t nkeep = -1;// Orthogonalize to at most the column dim 
    if (num > 0 && num <= ncols(M) && num <= nrows(M))
        {
        nkeep = num;
        }
    else
        {
        nkeep = std::min(nrows(M), ncols(M));
        }

    Vector dots(nkeep);
    MatrixRef Mcols;
    VectorRef dotsref, 
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

void
SVDRef(MatrixRefc const& M,
       MatrixRef  const& U, 
       VectorRef  const& D, 
       MatrixRef  const& V,
       Real thresh)
    {
    auto Mr = nrows(M), 
         Mc = ncols(M);

    if(Mr > Mc)
        {
        SVDRef(transpose(M),V,D,U,thresh);
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

    for(auto& el : D) // el = sqrt(fabs(el));
        {
        //This formula zeroes out any negative evals,
        //smaller error than using their abs value
        el = sqrt((fabs(el)+el)/2);
        }

    //Put result of Mt*U==(V*D) in V storage
    mult(transpose(M),U,V);
    for(auto c : index1(D)) 
        {
        if(D(c) > 0) column(V,c) /= D(c);
        }

    size_t start = 2;
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
    auto mv = move(rho);
    reduceCols(mv,n);
    mult(M,columns(V,start,Mr),mv);

    //b should be close to diagonal
    //but may not be perfect - fix it up below
    auto b = rows(transpose(U),start,Mr)*mv;
   
    auto d = subVector(D,start,Mr);
    Matrix u(n,n),
           v(n,n);
    SVDRef(b,u,d,v,thresh);

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

} //namespace itensor
