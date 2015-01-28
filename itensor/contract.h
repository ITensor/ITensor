//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_CONTRACT_H
#define __ITENSOR_CONTRACT_H

#include <array>

#include "print.h"
#include "permutation.h"
#include "simpletensor.h"

namespace itensor {

template<typename RangeT>
using RTref = tensorref<Real,RangeT>;
using Label = std::vector<int>;

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


template<typename R1, typename R2>
void 
reshape(const RTref<R1>& T, 
        const Permutation& P, 
        RTref<R2>& res);

template<typename RangeT>
void 
reshape(const RTref<RangeT>& T, 
        const Permutation& P, 
        tensor<Real,Range>& res);

//Callable is any function func(Real& x, Real y)
//default is func = [](Real& x, Real y) { x = y; };
template<typename R1, typename R2, typename Callable>
void 
reshape(const RTref<R1>& T, 
        const Permutation& P, 
        RTref<R2>& res,
        const Callable& func);

template<typename RangeT>
void 
contract(const RTref<RangeT>& A, const Label& ai, 
         const RTref<RangeT>& B, const Label& bi, 
         RTref<RangeT>& C,       const Label& ci);

template<typename RangeT>
void 
contractloop(const RTref<RangeT>& A, const Label& ai, 
             const RTref<RangeT>& B, const Label& bi, 
             RTref<RangeT>& C,       const Label& ci,
             const Args& args = Global::args());

template<typename T>
void
printv(const std::vector<T>& t)
    {
    print("{ ");
    for(const auto& i : t) print(i," ");
    println("}");
    }
template<typename T,typename F>
void
printv(const std::vector<T>& t,
       const F& f)
    {
    print("{ ");
    for(const auto& i : t) 
        {
        f(i);
        print(" ");
        }
    println("}");
    }
template<typename T, size_t size>
void
printv(const std::array<T,size>& t)
    {
    print("{ ");
    for(const auto& i : t) print(i," ");
    println("}");
    }
template<typename T>
void
printv(const autovector<T>& t)
    {
    print("{ ");
    for(auto i = t.mini(); i <= t.maxi(); ++i)
        {
        print(t.fast(i)," ");
        }
    println("}");
    }
#define PRI(a) print(#a,": "); printv(a);
#define PRIL(a,l) print(#a,": "); printv(a,l);

inline 
std::ostream& 
operator<<(std::ostream& s, const Label& A)
    {
    for(const auto& a : A) s << a << " ";
    s << "\n";
    return s;
    }

template<typename T>
long 
findIndex(const std::vector<T>& v, 
          const T& t)
    {
    using size_type = typename std::vector<T>::size_type;
    for(size_type i = 0; i < v.size(); ++i)
        if(v[i] == t) return i;
    return -1;
    }

template<typename R>
inline 
Real 
dist(const RTref<R>& A, const RTref<R>& Ach)
    {
    Real dif = 0.0;
    for(long i = 0; i < A.size(); i++)
        dif += sqr(A.v(i) - Ach.v(i));
    return std::sqrt(dif);
    };



///
/// Implementations
///

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
    long ncont = 1;
    for(long i = 0; i < rL; ++i)
    for(long j = 0; j < rR; ++j)
        if(Lis[i] == Ris[j])
            {
            //Negative entries in 
            //Lind, Rind indicate
            //contracted indices
            Lind[i] = -ncont;
            Rind[j] = -ncont;
            checkCont(Lis[i],Ris[j]);
            ++ncont;
            break;
            }

    //Go through and assign uncontracted entries of Lind,Rind
    //the integers ncont+1,ncont+2,...
    auto uu = ncont;
    for(long j = 0; j < rL; ++j)
        {
        if(Lind[j] == 0) Lind[j] = uu++;
        }
    for(long j = 0; j < rR; ++j)
        {
        if(Rind[j] == 0) Rind[j] = uu++;
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

template<typename R1, typename R2, typename Callable>
void 
reshape(const RTref<R1>& T, 
        const Permutation& P, 
        RTref<R2>& res,
        const Callable& func)
    {
    auto r = P.size();

#ifdef DEBUG
    if(res.size() != T.size()) Error("Mismatched storage sizes in reshape");
#endif

    //find largest dimension of T,
    //size "big" and position "bigind"
    long bigind = 0, 
         big = T.n(0);
    for(int j = 1; j < r; ++j)
        if(big < T.n(j))
            {
            big = T.n(j); 
            bigind = j;
            }

    auto stept = T.stride(bigind);
    auto stepr = res.stride(P.dest(bigind));

    detail::GCounter c(0,r-1,0);
    for(int i = 0; i < r; ++i)
        c.setInd(i,0,T.n(i)-1);
    c.setInd(bigind,0,0);		// This one we leave at 0 only

    Label Ti(r), 
          ri(r);
    for(; c.notDone(); ++c)
        {
        for(int i = 0; i < r; ++i)
            Ti[i] = ri[P.dest(i)] = c.i.fast(i);

        auto* pr = &res.vref(ind(res,ri));
        auto* pt = &T.v(ind(T,Ti));
        for(int k = 0; k < big; ++k)
            {
            func(*pr,*pt);
            //*pr = *pt;
            pr += stepr;
            pt += stept;
            }
        }
    }

namespace detail {
template<typename T>
void 
assign(T& r1, T r2) { r1 = r2; }
template<typename T>
void
plusEq(T& r1, T r2) { r1 += r2; }
};

template<typename R1, typename R2>
void 
reshape(const RTref<R1>& T, 
        const Permutation& P, 
        RTref<R2>& res)
    {
    reshape(T,P,res,detail::assign<Real>);
    }

template<typename RangeT>
void 
reshape(const RTref<RangeT>& T, 
        const Permutation& P, 
        tensor<Real,Range>& res)
    {
    auto r = P.size();
    std::vector<long> resdims(r);
    for(int i = 0; i < r; ++i)
        resdims[P.dest(i)] = T.n(i);
    res.resize(resdims);
    tensorref<Real,Range> res_ref(res.data(),res.inds());
    reshape(T,P,res_ref);
    }

}; //namespace itensor

#endif
