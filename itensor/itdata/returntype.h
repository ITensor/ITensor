//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
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
    using Test = std::invoke_result_t<TestRet<Task,frontType<TList>>>;
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

