//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_CONTRACT_H
#define __ITENSOR_CONTRACT_H

#include "permute.h"
#include "global.h"

namespace itensor {

template<typename Inds, typename Func>
long
computeLabels(const Inds& Lis,
              long rL,
              const Inds& Ris,
              long rR,
              Label& Lind,
              Label& Rind,
              const Func& checkCont);

template<typename Inds>
long
computeLabels(const Inds& Lis,
              long rL,
              const Inds& Ris,
              long rR,
              Label& Lind,
              Label& Rind);


template<typename RangeT>
void 
contract(TenRefc<RangeT> A, const Label& ai, 
         TenRefc<RangeT> B, const Label& bi, 
         TenRef<RangeT>  C, const Label& ci);

template<typename RangeT>
void 
contract(const Tensor<RangeT>& A, const Label& ai, 
         const Tensor<RangeT>& B, const Label& bi, 
         Tensor<RangeT>&  C, const Label& ci);

template<typename RangeT>
void 
contractloop(TenRefc<RangeT> A, const Label& ai, 
             TenRefc<RangeT> B, const Label& bi, 
             TenRef<RangeT>  C, const Label& ci,
             const Args& args = Global::args());

template<typename RangeT>
void 
contractloop(const Tensor<RangeT>& A, const Label& ai, 
             const Tensor<RangeT>& B, const Label& bi, 
             Tensor<RangeT>&  C, const Label& ci,
             const Args& args = Global::args());


///
/// Implementations
///

template<typename RangeT>
void 
contract(const Tensor<RangeT>& A, const Label& ai, 
         const Tensor<RangeT>& B, const Label& bi, 
         Tensor<RangeT>&  C, const Label& ci)
    {
    contract(makeRef(A),ai,makeRef(B),bi,makeRef(C),ci);
    }

template<typename RangeT>
void 
contractloop(const Tensor<RangeT>& A, const Label& ai, 
             const Tensor<RangeT>& B, const Label& bi, 
             Tensor<RangeT>&  C, const Label& ci,
             const Args& args)
    {
    contractloop(makeRefc(A),ai,makeRefc(B),bi,makeRef(C),ci,args);
    }

template<typename Inds, typename Func>
long
computeLabels(const Inds& Lis,
              long rL,
              const Inds& Ris,
              long rR,
              Label& Lind,
              Label& Rind,
              const Func& checkCont)
    {
    //Set Lind, Rind to zero. Special value 0 marks
    //uncontracted indices. Later will assign unique numbers
    //to these entries in Lind and Rind
    Lind.assign(rL,0);
    Rind.assign(rR,0);

    //Count number of contracted indices,
    //set corresponding entries of Lind, Rind
    //to 1,2,...,ncont
       // if(Lis[i] == Ris[j])
#ifdef DEBUG
    if(Global::debug2()) print("Contracting: ");
#endif
    long ncont = 0;
    for(long i = 0; i < rL; ++i)
    for(long j = 0; j < rR; ++j)
        if(Lis[i] == Ris[j])
            {
#ifdef DEBUG
            if(Global::debug2())
                {
                print(Lis[i],", ");
                }
#endif
            //Negative entries in 
            //Lind, Rind indicate
            //contracted indices
            Lind[i] = -(1+ncont);
            Rind[j] = -(1+ncont);
            checkCont(Lis[i],Ris[j]);
            ++ncont;
            break;
            }
#ifdef DEBUG
    if(Global::debug2()) 
        {
        println(" p = ",Lis.r()+Ris.r()-ncont);
        }
#endif

    //Go through and assign uncontracted entries of Lind,Rind
    //the integers ncont+1,ncont+2,...
    auto uu = ncont;
    for(long j = 0; j < rL; ++j)
        {
        if(Lind[j] == 0) Lind[j] = ++uu;
        }
    for(long j = 0; j < rR; ++j)
        {
        if(Rind[j] == 0) Rind[j] = ++uu;
        }
    return ncont;
    }

template<typename Inds>
long
computeLabels(const Inds& Lis,
              long rL,
              const Inds& Ris,
              long rR,
              Label& Lind,
              Label& Rind)
    {
    using ind = typename Inds::value_type;
    auto nocheck = [](const ind& li,const ind& ri) { };
    return computeLabels(Lis,rL,Ris,rR,Lind,Rind,nocheck);
    }

} //namespace itensor

#endif
