//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_RETURNTYPE_H
#define __ITENSOR_RETURNTYPE_H

#include "itensor/util/stdx.h"
#include "itensor/util/typelist.h"
#include "itensor/itdata/itdata.h"

namespace itensor {

namespace detail { 

//Converts (Ret==void)->NoneType, otherwise keeps Ret->Ret
template<typename Ret>
using RetOrNone = stdx::conditional_t<std::is_void<Ret>::value,NoneType,Ret>;

template<typename Task, typename Storage>
NoneType
testRetImpl(stdx::choice<3>, Task& t, Storage& s, ManageStore& m)
    {
    return NoneType{};
    }
template<typename Task, typename Storage>
auto 
testRetImpl(stdx::choice<2>, Task& t, Storage& s, ManageStore& m)
    -> RetOrNone<decltype(doTask(t,s))>
    {
    using Ret = RetOrNone<decltype(doTask(t,s))>;
    return Ret{};
    }
template<typename Task, typename Storage>
auto 
testRetImpl(stdx::choice<1>, Task& t, Storage& s, ManageStore& m)
    -> RetOrNone<decltype(doTask(t,s,m))>
    {
    using Ret = RetOrNone<decltype(doTask(t,s,m))>;
    return Ret{};
    }

//template<typename Task, typename Storage>
//struct DoTaskTester
//    {
//    auto
//    operator()(Task & t, Storage & s, ManageStore & m)
//        {
//        return doTask(t,s,m);
//        }
//    };

template<typename Task, typename Storage>
struct TestRet
    {
    Task& t;
    Storage& s;
    ManageStore& m;
    TestRet(Task& t_, Storage& s_, ManageStore& m_)
        : t(t_), s(s_), m(m_) 
        { }
    auto
    operator()() -> decltype(testRetImpl(stdx::select_overload{},t,s,m))
        {
        return testRetImpl(stdx::select_overload{},t,s,m);
        }
    };

//GetRType walks the typelist TList of storage types until
//TestRet finds a non-trivial implementation of doTask
template<typename Task, typename TList>
struct GetRType : GetRType<Task,popFront<TList>>
    {
    using Test = stdx::result_of_t<TestRet<Task,frontType<TList>>()>;
    using Parent = GetRType<Task,popFront<TList>>;
    using RType = stdx::conditional_t<not std::is_same<Test,NoneType>::value,
                                     Test,
                                     typename Parent::RType>;
    };
template<typename Task>
struct GetRType<Task,TypeList<>>
    {
    using RType = NoneType;
    };

} //namespace detail

template<typename Task, typename TList>
using DoTaskReturn = typename detail::GetRType<stdx::remove_reference_t<Task>,TList>::RType;

} //namespace itensor


#endif

