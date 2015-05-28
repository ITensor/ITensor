//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_DOTASK_TEMPLATES_H
#define __ITENSOR_DOTASK_TEMPLATES_H

#include "itensor/itdata/itdata.h"
#include "itensor/detail/call_rewrite.h"

namespace itensor {

template<typename F>
void
doTask(GenerateIT<F>& G, const ITReal& d, ManagePtr& mp)
    { 
    if(G.isComplex())
        {
        auto* nd = mp.makeNewData<ITCplx>(d.size());
        for(auto j = 0ul; j < nd->csize(); ++j)
            nd->set(j,G.f());
        }
    else
        {
        auto* pd = mp.modifyData(d);
        std::generate(pd->begin(),pd->end(),[&G](){ return std::real(G.f()); });
        }
    }

template<typename F>
void
doTask(GenerateIT<F>& G, const ITCplx& d, ManagePtr& mp)
    { 
    if(G.isComplex())
        {
        auto* pd = mp.modifyData(d);
        for(auto j = 0ul; j < pd->csize(); ++j)
            pd->set(j,G.f());
        }
    else
        {
        auto* nd = mp.makeNewData<ITReal>(d.csize());
        std::generate(nd->begin(),nd->end(),[&G](){ return std::real(G.f()); });
        }
    }

template<typename F>
void
doTask(GenerateIT<F>& G, const IQTData& cd, ManagePtr& mp)
    {
    if(G.isComplex())
        {
        Error("Complex version of IQTensor generate not yet supported");
        }
    else
        {
        auto* pd = mp.modifyData(cd);
        std::generate(pd->data.begin(),pd->data.end(),[&G](){return std::real(G.f()); });
        }
    }

template<typename F>
void
doTask(ApplyIT<F>& A, ITReal& d)
    { 
    for(auto& elt : d) elt = detail::call<Real>(A.f,elt);
    }

template<typename F>
void
doTask(ApplyIT<F>& A, ITCplx& d)
    { 
    for(auto j = 0ul; j < d.csize(); ++j)
        {
        auto res = detail::call<Cplx>(A.f,d.get(j));
        d.set(j,res);
        }
    }

template <typename F, typename T,
          typename std::enable_if<std::is_same<T,std::result_of_t<F(T)>>::value>::type* = nullptr>
void
doTask(ApplyIT<F>& A, ITDiag<T>& d) 
    { 
    if(d.allSame()) 
        {
        d.val = detail::call<T>(A.f,d.val);
        }
    else
        {
        for(auto& elt : d.store) elt = detail::call<T>(A.f,elt);
        }
    }


template <typename F>
void
doTask(ApplyIT<F>& A, IQTData& d)
    {
    for(auto& elt : d.data)
        elt = A.f(elt);
    }

template<typename F>
void
doTask(VisitIT<F>& V, const ITReal& d)
    { 
    for(auto& elt : d) detail::call<void>(V.f,V.scale_fac * elt);
    }

template<typename F>
void
doTask(VisitIT<F>& V, const ITCplx& d)
    { 
    for(auto j = 0ul; j < d.csize(); ++j)
        {
        detail::call<void>(V.f,V.scale_fac * d.get(j));
        }
    }

template <typename F>
void
doTask(VisitIT<F>& V, const IQTData& d)
    {
    for(const auto& elt : d.data)
        V.f(elt*V.scale_fac);
    }

} // namespace itensor

#endif
