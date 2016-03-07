//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_APPLYFUNC_H
#define __ITENSOR_APPLYFUNC_H

#include "itensor/itdata/dotask.h"

namespace itensor {

namespace detail {

template<typename FuncObj>
struct FuncHolder
    {
    FuncObj& f;
    FuncHolder(FuncObj&& f_) : f(f_) { }
    };

template<typename F, typename Storage>
void
applyFunc_impl(stdx::choice<5>, FuncHolder<F>& H, Storage const& s, ManageStore& m)
    {
    throw ITError("applyFunc: function object has no operator() method for storage type");
    }
template<typename F, typename Storage>
auto
applyFunc_impl(stdx::choice<4>, FuncHolder<F>& H, Storage const& s, ManageStore& m)
    -> stdx::enable_if_t<not std::is_void<decltype(H.f(*m.modifyData(s)))>::value,decltype((H.f(*m.modifyData(s))))>
    {
    auto *ncs = m.modifyData(s);
    return H.f(*ncs);
    }
template<typename F, typename Storage>
auto
applyFunc_impl(stdx::choice<3>, FuncHolder<F>& H, Storage const& s, ManageStore& m)
    -> stdx::enable_if_t<std::is_void<decltype(H.f(*m.modifyData(s)))>::value,void>
    {
    auto *ncs = m.modifyData(s);
    H.f(*ncs);
    }
template<typename F, typename Storage>
auto
applyFunc_impl(stdx::choice<2>, FuncHolder<F>& H, Storage const& s, ManageStore& m) 
    -> stdx::enable_if_t<not std::is_void<decltype(H.f(s))>::value,decltype(H.f(s))>
    {
    return H.f(s);
    }
template<typename F, typename Storage>
auto
applyFunc_impl(stdx::choice<1>, FuncHolder<F>& H, Storage const& s, ManageStore& m) 
    -> stdx::enable_if_t<std::is_void<decltype(H.f(s))>::value,void>
    {
    H.f(s);
    }

} //namespace detail


template<typename FuncObj, typename Storage>
auto
doTask(detail::FuncHolder<FuncObj> & H, Storage const& s, ManageStore & m) 
    -> decltype(detail::applyFunc_impl(stdx::select_overload{},H,s,m))
    { 
    return detail::applyFunc_impl(stdx::select_overload{},H,s,m);
    }

template<typename F>
auto
applyFunc(F&& f, PData & store)
    -> decltype(doTask(detail::FuncHolder<F>{std::forward<F>(f)},store))
    {
    return doTask(detail::FuncHolder<F>{std::forward<F>(f)},store);
    }

template<typename F>
auto
applyFunc(F&& f, CPData const& store)
    -> decltype(doTask(detail::FuncHolder<F>{std::forward<F>(f)},store))
    {
    return doTask(detail::FuncHolder<F>{std::forward<F>(f)},store);
    }

} //namespace itensor

#endif
