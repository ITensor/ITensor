//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#pragma once

#include <array>

#include "print.h"
#include "permutation.h"
#include "simpletensor.h"

namespace itensor {

template<typename RangeT>
using RTref = tensorref<Real,RangeT>;
using Label = std::vector<int>;

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
    for(const auto& i : t) print(i," ");
    println();
    }
template<typename T, size_t size>
void
printv(const std::array<T,size>& t)
    {
    for(const auto& i : t) print(i," ");
    println();
    }
template<typename T>
void
printv(const autovector<T>& t)
    {
    for(auto i = t.mini(); i <= t.maxi(); ++i)
        {
        print(t.fast(i)," ");
        }
    println();
    }
#define PRI(a) print(#a,": "); printv(a);

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
void inline
assign(Real& r1, Real r2) { r1 = r2; }
};

template<typename R1, typename R2>
void 
reshape(const RTref<R1>& T, 
        const Permutation& P, 
        RTref<R2>& res)
    {
    reshape(T,P,res,detail::assign);
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
