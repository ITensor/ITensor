//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_CONTRACT_H
#define __ITENSOR_CONTRACT_H

#include "itensor/tensor/vec.h"
#include "itensor/util/args.h"
#include "itensor/util/range.h"
#include "itensor/detail/gcounter.h"

namespace itensor {

template<typename RangeT, typename VA, typename VB>
void 
contract(TenRefc<RangeT,VA> A, Labels const& ai, 
         TenRefc<RangeT,VB> B, Labels const& bi, 
         TenRef<RangeT,common_type<VA,VB>>  C, 
         Labels const& ci,
         Real alpha = 1.,
         Real beta = 0.);

template<typename R, typename VA, typename VB>
void 
contract(Ten<R,VA> const& A, Labels const& ai, 
         Ten<R,VB> const& B, Labels const& bi, 
         Ten<R,common_type<VA,VB>> & C, Labels const& ci,
         Real alpha = 1.,
         Real beta = 0.);

template<typename range_type>
void 
contractloop(TenRefc<range_type> A, Labels const& ai, 
             TenRefc<range_type> B, Labels const& bi, 
             TenRef<range_type>  C, Labels const& ci,
             Args const& args = Args::global());

template<typename range_type>
void 
contractloop(Ten<range_type> const& A, Labels const& ai, 
             Ten<range_type> const& B, Labels const& bi, 
             Ten<range_type>      & C, Labels const& ci,
             Args const& args = Args::global());


//All indices of B contracted
//(A can have some uncontracted indices)
template<typename DiagElsA, typename RangeT, typename VB, typename VC>
void 
contractDiagFull(DiagElsA           const& A, Labels const& ai, 
                 TenRefc<RangeT,VB> const& B, Labels const& bi, 
                 VecRef<VC>         const& C, Labels const& ci,
                 IntArray astarts = IntArray{});

//Some indices of B uncontracted
template<typename DiagElsA, typename RangeT, typename VB, typename VC>
void 
contractDiagPartial(DiagElsA           const& A, Labels const& ai,
                    TenRefc<RangeT,VB> const& B, Labels const& bi, 
                    TenRef<RangeT,VC>  const& C, Labels const& ci,
                    IntArray                  astarts = IntArray{});

//Non-contracting product
template<class TA, class TB, class TC>
void 
ncprod(TA && A, Labels const& ai, 
       TB && B, Labels const& bi, 
       TC && C, Labels const& ci);

template<typename Inds, typename Func>
long
computeLabels(Inds const& Lis,
              long rL,
              Inds const& Ris,
              long rR,
              Labels& Lind,
              Labels& Rind,
              Func const& checkCont);

template<typename Inds>
long
computeLabels(Inds const& Lis,
              long rL,
              Inds const& Ris,
              long rR,
              Labels & Lind,
              Labels & Rind);

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

template<typename R, typename VA, typename VB>
void 
contract(Ten<R,VA> const& A, Labels const& ai, 
         Ten<R,VB> const& B, Labels const& bi, 
         Ten<R,common_type<VA,VB>> & C, Labels const& ci,
         Real alpha,
         Real beta)
    {
    contract(makeRef(A),ai,makeRef(B),bi,makeRef(C),ci,alpha,beta);
    }

template<typename range_type>
void
contractloop(Ten<range_type> const& A, Labels const& ai, 
             Ten<range_type> const& B, Labels const& bi, 
             Ten<range_type>      & C, Labels const& ci,
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
              Labels& Lind,
              Labels& Rind,
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
              Labels & Lind,
              Labels & Rind)
    {
    using ind = typename Inds::value_type;
    auto nocheck = [](const ind& li,const ind& ri) { };
    return computeLabels(Lis,rL,Ris,rR,Lind,Rind,nocheck);
    }

//contractDiagPartial:
//Some indices of B uncontracted
//DiagElsA is any function object returning
//diagonal elements of A (such as a VectorRefc)
//
//Although all indices of A proceed in lockstep,
//some may not start at zero, but at various values
//provided by "astarts"; this is useful for blocks
//of diagonal IQTensors (QDiag storage) where the
//diagonal elements do not necessarily start in 
//the 0,0,0,...,0 entry of a particular block
template<typename DiagElsA, typename RangeT, typename VB, typename VC>
void 
contractDiagPartial(DiagElsA           const& A, Labels const& al,
                    TenRefc<RangeT,VB> const& B, Labels const& bl, 
                    TenRef<RangeT,VC>  const& C, Labels const& cl,
                    IntArray astarts)
    {
    if(astarts.empty()) astarts.assign(al.size(),0);
    auto bstart = 0ul;
    auto cstart = 0ul;
    auto b_cstride = 0ul; //B contracted stride
    int nbu = 0;          //# B uncont. inds.
    for(auto ib : range(bl))
        {
        auto ia = find_index(al,bl[ib]);
        if(ia >= 0)
            {
            b_cstride += B.stride(ib);
            bstart += astarts[ia]*B.stride(ib);
            }
        else       
            {
            nbu += 1;
            }
        }

    auto c_cstride = 0ul;
    for(auto ic : range(cl))
        {
        auto ia = find_index(al,cl[ic]);
        if(ia >= 0) 
            {
            c_cstride += C.stride(ic);
            cstart += astarts[ia]*C.stride(ic);
            }
        }

    auto bustride = IntArray(nbu,0);
    auto custride = IntArray(nbu,0);
    auto GC = detail::GCounter(nbu);
    int n = 0;
    for(auto ib : range(bl))
        {
        if(bl[ib] > 0) //uncontracted
            {
#ifdef DEBUG
            if(n >= nbu) Error("n out of range");
#endif
            GC.setRange(n,0,B.extent(ib)-1);
            bustride[n] = B.stride(ib);
            auto ic = find_index(cl,bl[ib]);
#ifdef DEBUG
            if(ic < 0) Error("Index not found");
#endif
            custride[n] = C.stride(ic);
            ++n;
            }
        }
    auto pb = MAKE_SAFE_PTR(B.data(),B.size());
    auto pc = MAKE_SAFE_PTR(C.data(),C.size());
    for(;GC.notDone();++GC)
        {
        size_t coffset = 0;
        size_t boffset = 0;
        for(auto i : range(nbu))
            {
            auto ii = GC[i];
            boffset += ii*bustride[i];
            coffset += ii*custride[i];
            }
        for(auto J : range(A))
            {
            pc[cstart+J*c_cstride+coffset] += A(J)*pb[bstart+J*b_cstride+boffset];
            }
        }
    }

// C = A*B
//case where all of B's indices are contracted with A,
//making C diagonal
template<typename DiagElsA, typename RangeT, typename VB, typename VC>
void 
contractDiagFull(DiagElsA           const& A, Labels const& al, 
                 TenRefc<RangeT,VB> const& B, Labels const& bl, 
                 VecRef<VC>         const& C, Labels const& cl,
                 IntArray astarts)
    {
    if(astarts.empty()) astarts.assign(al.size(),0);
    auto bstart = 0ul;
    long b_cstride = 0; //total stride of contracted inds of B
    for(auto ib : range(bl))
        {
        auto ia = find_index(al,bl[ib]);
        if(ia >= 0) 
            {
            b_cstride += B.stride(ib);
            bstart += astarts[ia]*B.stride(ib);
            }
        }

    auto pb = MAKE_SAFE_PTR(B.data(),B.size());
    if(C.size() == 1)
        {
        auto *Cval = C.data();
        for(auto J : range(A))
            {
            *Cval += A(J)*pb[bstart+J*b_cstride];
            }
        }
    else
        {
        auto pc = MAKE_SAFE_PTR(C.data(),C.size());
        for(auto J : range(A))
            {
            pc[J] += A(J)*pb[bstart+J*b_cstride];
            }
        }
    }

} //namespace itensor

#include "itensor/tensor/contract_impl.h"

#endif
