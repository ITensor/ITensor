//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_CONTRACT_H
#define __ITENSOR_CONTRACT_H

#include "itensor/tensor/vec.h"
#include "itensor/util/args.h"
#include "itensor/detail/gcounter.h"

namespace itensor {

template<typename RangeT>
void 
contract(TenRefc<RangeT> A, Label const& ai, 
         TenRefc<RangeT> B, Label const& bi, 
         TenRef<RangeT>  C, Label const& ci);

void 
contract(Tensor const& A, Label const& ai, 
         Tensor const& B, Label const& bi, 
         Tensor      & C, Label const& ci);

template<typename range_type>
void 
contractloop(TenRefc<range_type> A, Label const& ai, 
             TenRefc<range_type> B, Label const& bi, 
             TenRef<range_type>  C, Label const& ci,
             Args const& args = Args::global());

template<typename range_type>
void 
contractloop(Ten<range_type> const& A, Label const& ai, 
             Ten<range_type> const& B, Label const& bi, 
             Ten<range_type>      & C, Label const& ci,
             Args const& args = Args::global());


//All indices of B contracted
//(A can have some uncontracted indices)
template<typename DiagElsA, typename RangeT>
void 
contractDiagFull(DiagElsA const& A, Label const& ai, 
                 TenRefc<RangeT> B, Label const& bi, 
                 VectorRef          C, Label const& ci);

//Some indices of B uncontracted
template<typename DiagElsA, typename RangeT>
void 
contractDiagPartial(DiagElsA const& A,  Label const& ai,
                    TenRefc<RangeT> B, Label const& bi, 
                    TenRef<RangeT>  C, Label const& ci);

template<typename Inds, typename Func>
long
computeLabels(Inds const& Lis,
              long rL,
              Inds const& Ris,
              long rR,
              Label& Lind,
              Label& Rind,
              Func const& checkCont);

template<typename Inds>
long
computeLabels(Inds const& Lis,
              long rL,
              Inds const& Ris,
              long rR,
              Label & Lind,
              Label & Rind);

template<typename T>
long 
find_index(std::vector<T> const& v, 
           T const& t)
    {
    for(decltype(v.size()) i = 0; i < v.size(); ++i)
        if(v[i] == t) return i;
    return -1;
    }
template<typename T, size_t MaxSize>
long 
find_index(VarArray<T,MaxSize> const& v, 
           T const& t)
    {
    for(decltype(v.size()) i = 0; i < v.size(); ++i)
        if(v[i] == t) return i;
    return -1;
    }
template<typename T, size_t MaxSize>
long 
find_index(InfArray<T,MaxSize> const& v, 
           T const& t)
    {
    for(decltype(v.size()) i = 0; i < v.size(); ++i)
        if(v[i] == t) return i;
    return -1;
    }


///
/// Implementations
///

void inline
contract(Tensor const& A, Label const& ai, 
         Tensor const& B, Label const& bi, 
         Tensor      & C, Label const& ci)
    {
    contract(makeRef(A),ai,makeRef(B),bi,makeRef(C),ci);
    }

template<typename range_type>
void
contractloop(Ten<range_type> const& A, Label const& ai, 
             Ten<range_type> const& B, Label const& bi, 
             Ten<range_type>      & C, Label const& ci,
             Args const& args)
    {
    contractloop(makeRefc(A),ai,makeRefc(B),bi,makeRef(C),ci,args);
    }

template<typename Inds, typename Func>
long
computeLabels(Inds const& Lis,
              long rL,
              Inds const& Ris,
              long rR,
              Label& Lind,
              Label& Rind,
              Func const& checkCont)
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
    long ncont = 0;
    for(long i = 0; i < rL; ++i)
    for(long j = 0; j < rR; ++j)
        if(Lis[i] == Ris[j])
            {
            //Negative entries in 
            //Lind, Rind indicate
            //contracted indices
            Lind[i] = -(1+ncont);
            Rind[j] = -(1+ncont);
            checkCont(Lis[i],Ris[j]);
            ++ncont;
            break;
            }

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
computeLabels(Inds const& Lis,
              long rL,
              Inds const& Ris,
              long rR,
              Label & Lind,
              Label & Rind)
    {
    using ind = typename Inds::value_type;
    auto nocheck = [](const ind& li,const ind& ri) { };
    return computeLabels(Lis,rL,Ris,rR,Lind,Rind,nocheck);
    }

//Some indices of B uncontracted
//DiagElsA is any function object returning
//diagonal elements of A (such as a VectorRefc)
template<typename DiagElsA, typename RangeT>
void 
contractDiagPartial(DiagElsA const& A,  Label const& ai,
                    TenRefc<RangeT> B, Label const& bi, 
                    TenRef<RangeT>  C, Label const& ci)
    {
    using A_size_type = decltype(A.size());
    size_t b_cstride = 0; //B contracted stride
    int nbu = 0;          //# B uncont. inds.
    for(auto j = 0ul; j < bi.size(); ++j)
        {
        //if index j is contracted, add its stride to n_cstride:
        if(bi[j] < 0) b_cstride += B.stride(j);
        else            ++nbu;
        }

    long a_ustride = 0; //total stride of uncontracted
                        //inds of A (infer from C)
    for(auto i = 0ul; i < ci.size(); ++i)
        {
        auto j = find_index(ai,ci[i]);
        if(j >= 0) a_ustride += C.stride(i);
        }

    Label bstride(nbu,0),
          cstride(nbu,0);
    detail::GCounter GC(nbu);
    int n = 0;
    for(decltype(bi.size()) j = 0; j < bi.size(); ++j)
        {
        if(bi[j] > 0)
            {
#ifdef DEBUG
            if(n >= nbu) Error("n out of range");
#endif
            GC.setRange(n,0,B.extent(j)-1);
            bstride[n] = B.stride(j);
            auto k = find_index(ci,bi[j]);
#ifdef DEBUG
            if(k < 0) Error("Index not found");
#endif
            cstride[n] = C.stride(k);
            ++n;
            }
        }
    auto pb = MAKE_SAFE_PTR(B.data(),B.size());
    auto pc = MAKE_SAFE_PTR(C.data(),C.size());
    for(;GC.notDone();++GC)
        {
        size_t coffset = 0,
               boffset = 0;
        for(auto i = 0; i < nbu; ++i)
            {
            auto ii = GC[i];
            boffset += ii*bstride[i];
            coffset += ii*cstride[i];
            }
        for(A_size_type J = 0; J < A.size(); ++J)
            {
            pc[J*a_ustride+coffset] += A(1+J)*pb[J*b_cstride+boffset];
            }
        }
    }

template<typename DiagElsA, typename RangeT>
void 
contractDiagFull(DiagElsA const& A,  Label const& ai, 
                 TenRefc<RangeT> B, Label const& bi, 
                 VectorRef          C,  Label const& ci)
    {
    using A_size_type = decltype(A.size());

    long b_cstride = 0; //total stride of contracted inds of B
    for(auto j = 0ul; j < bi.size(); ++j)
        if(bi[j] < 0) b_cstride += B.stride(j);

    auto pb = MAKE_SAFE_PTR(B.data(),B.size());
    if(C.size() == 1)
        {
        auto *Cval = C.data();
        for(A_size_type J = 0; J < A.size(); ++J)
            *Cval += A(1+J)*pb[J*b_cstride];
        }
    else
        {
        auto pc = MAKE_SAFE_PTR(C.data(),C.size());
        for(A_size_type J = 0; J < A.size(); ++J)
            pc[J] += A(1+J)*pb[J*b_cstride];
        }
    }

} //namespace itensor

#endif
