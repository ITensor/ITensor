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
#ifndef __ITENSOR_TYPELIST_H
#define __ITENSOR_TYPELIST_H

//
// Ideas for improvement:
// o Add indexOf (returns index of type T if containsType<TList,T>==true,
//   such that getType<TList,indexOf<TList,T>()> == T)
// o Add insert
//

namespace itensor {

template<typename... Ts>
struct TypeList 
    { 
    size_t static constexpr
    size() { return 0ul; }

    template<typename... NewTypes>
    using pushFront = TypeList<NewTypes...,Ts...>;

    template<typename... NewTypes>
    using pushBack = TypeList<Ts...,NewTypes...>;
    };

template<typename T, typename... Ts>
struct TypeList<T,Ts...> : TypeList<Ts...>
    {
    using Type = T;
    using Next = TypeList<Ts...>;

    size_t static constexpr
    size() { return 1ul+Next::size(); }

    template<typename... NewTypes>
    using pushFront = TypeList<NewTypes...,T,Ts...>;

    template<typename... NewTypes>
    using pushBack = TypeList<T,Ts...,NewTypes...>;
    };

//
// pushFront
//
template<typename TList, typename... NewTypes>
using pushFront = typename TList::template pushFront<NewTypes...>;

//
// pushBack
//
template<typename TList, typename... NewTypes>
using pushBack = typename TList::template pushBack<NewTypes...>;

//
// frontType
//
template<typename TList>
using frontType = typename TList::Type;

//
// popFront
//
template<typename TList>
using popFront = typename TList::Next;

struct NoneType { };

template<typename T, typename ElseType>
using ifTypeElse = typename std::conditional<not std::is_same<T,NoneType>::value,
                                      T,
                                      ElseType>::type;

template<typename T, typename TList>
struct CheckContainsType : CheckContainsType<T,popFront<TList>>
    {
    using Test = typename std::conditional<std::is_same<frontType<TList>,T>::value,T,NoneType>::type;
    using ParentResult = typename CheckContainsType<T,popFront<TList>>::Result;
    using Result = ifTypeElse<Test,ParentResult>;
    };
template<typename T>
struct CheckContainsType<T,TypeList<>>
    {
    using Result = NoneType;
    };

//
// containsType
//
template<typename TList, typename T, typename C = typename CheckContainsType<T,TList>::Result>
struct containsType : std::true_type 
    {
    constexpr operator bool() const noexcept { return true; }
    };

template<typename TList, typename T>
struct containsType<TList,T,NoneType> : std::false_type 
    { 
    constexpr operator bool() const noexcept { return false; }
    };


// May be nicer to use variadic template-template param
// see: http://ow.ly/OIurr

//template<size_t s, typename... Ts>
//struct PushFront { };
//
//template<typename TL, typename N, typename... Ts>
//struct PushFront<0ul,TL,N,Ts...>
//    {
//    using Result = TypeList<N,Ts...>;
//    };
//
//template<size_t s, typename TL, typename... Ts>
//struct PushFront<s,TL,Ts...> : PushFront<popFront<TL>::size(),popFront<TL>,Ts...,frontType<TL>>
//    {
//    using Next = PushFront<popFront<TL>::size(),popFront<TL>,Ts...,frontType<TL>>;
//    using Result = typename Next::Result;
//    };
//
////
//// pushFront
////
//template<typename TList, typename... NewTypes>
//using pushFront = typename PushFront<TList::size(),TList,NewTypes...>::Result;



namespace detail {
template<typename TL, size_t n, bool ok>
struct GetType : GetType<popFront<TL>,n-1,ok>
    {
    using Next = GetType<popFront<TL>,n-1,ok>;
    using Result = typename Next::Result;
    };
template<typename TL, bool ok>
struct GetType<TL,0ul,ok>
    {
    using Result = frontType<TL>;
    };
template<typename TL, size_t n>
struct GetType<TL,n,false>
    {
    using Result = NoneType;
    static_assert(n < TL::size(),"getType argument out of range");
    };
} //namespace detail

//
// getType
//
template<typename TList, size_t n>
using getType = typename detail::GetType<TList,n,n < TList::size()>::Result;

} //namespace itensor

#endif
