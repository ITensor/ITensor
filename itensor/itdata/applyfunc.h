//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_APPLYFUNC_H
#define __ITENSOR_APPLYFUNC_H

#include "itensor/itdata/dotask.h"

namespace itensor {

namespace detail {

template<typename F,typename Ret>
struct ApplyFunc
    { 
    using function_type = F;
    F& f;
    Ret r;
    ApplyFunc(F&& f_) : f(f_) { }
    template<typename S>
    void
    operator()(S& s) 
        { 
        r = f(s); 
        }
    operator Ret() const { return std::move(r); }
    };

template<typename F>
struct ApplyFunc<F,void>
    { 
    using function_type = F;
    F& f;
    ApplyFunc(F&& f_) : f(f_) { }
    template<typename S>
    void
    operator()(S& s) { f(s); }
    };

template<typename F, typename R, typename Storage>
void
applyFunc_impl(stdx::choice<3>, ApplyFunc<F,R>& A, const Storage& s, ManageStore& m)
    {
    throw ITError("applyFunc: function object has no operator() method for storage type");
    }

template<typename F, typename R, typename Storage>
auto
applyFunc_impl(stdx::choice<2>, ApplyFunc<F,R>& A, const Storage& s, ManageStore& m)
    -> decltype(A.f(s), void())
    {
    Storage& ncs = m.modifyData();
    A(ncs);
    }

template<typename F, typename R, typename Storage>
auto
applyFunc_impl(stdx::choice<1>, ApplyFunc<F,R>& A, const Storage& s, ManageStore& m) 
    -> decltype(A.f(s), void())
    {
    A(s);
    }

} //namespace detail

template<typename F, typename R, typename Storage>
void
doTask(detail::ApplyFunc<F,R> & A, Storage const& s, ManageStore & m) 
    { 
    detail::applyFunc_impl(stdx::select_overload{},A,s,m);
    }


//template<typename F>
//F
//applyFunc(F&& f, PData & store)
//    {
//    doTask(detail::ApplyFunc<F,void>{std::forward<F>(f)},store);
//    return f;
//    }
//
//template<typename F>
//F
//applyFunc(F&& f, CPData const& store)
//    {
//    doTask(detail::ApplyFunc<F,void>{std::forward<F>(f)},store);
//    return f;
//    }

//template<typename F>
//auto
//applyFunc(F&& f, PData & store)
//    -> decltype(doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store))
//    {
//    return doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store);
//    }
//
//template<typename F>
//auto
//applyFunc(F&& f, CPData const& store)
//    -> decltype(doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store))
//    {
//    return doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store);
//    }


//template<typename F>
//F
//applyFunc(F&& f, PData & store)
//    {
//    doTask(detail::ApplyFunc<F,void>{std::forward<F>(f)},store);
//    return f;
//    }
//
//template<typename F>
//F
//applyFunc(F&& f, CPData const& store)
//    {
//    doTask(detail::ApplyFunc<F,void>{std::forward<F>(f)},store);
//    return f;
//    }

//template<typename F>
//auto
//applyFunc(F&& f, PData & store)
//    -> decltype(doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store))
//    {
//    return doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store);
//    }
//
//template<typename F>
//auto
//applyFunc(F&& f, CPData const& store)
//    -> decltype(doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store))
//    {
//    return doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store);
//    }
//template<typename F>
//F
//applyFunc(F&& f, PData & store)
//    {
//    doTask(detail::ApplyFunc<F,void>{std::forward<F>(f)},store);
//    return f;
//    }
//
//template<typename F>
//F
//applyFunc(F&& f, CPData const& store)
//    {
//    doTask(detail::ApplyFunc<F,void>{std::forward<F>(f)},store);
//    return f;
//    }

//template<typename F>
//auto
//applyFunc(F&& f, PData & store)
//    -> decltype(doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store))
//    {
//    return doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store);
//    }
//
//template<typename F>
//auto
//applyFunc(F&& f, CPData const& store)
//    -> decltype(doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store))
//    {
//    return doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store);
//    }

//template<typename F>
//F
//applyFunc(F&& f, PData & store)
//    {
//    doTask(detail::ApplyFunc<F,void>{std::forward<F>(f)},store);
//    return f;
//    }
//
//template<typename F>
//F
//applyFunc(F&& f, CPData const& store)
//    {
//    doTask(detail::ApplyFunc<F,void>{std::forward<F>(f)},store);
//    return f;
//    }

//template<typename F>
//auto
//applyFunc(F&& f, PData & store)
//    -> decltype(doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store))
//    {
//    return doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store);
//    }
//
//template<typename F>
//auto
//applyFunc(F&& f, CPData const& store)
//    -> decltype(doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store))
//    {
//    return doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store);
//    }

} //namespace itensor

#endif
