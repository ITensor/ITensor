//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TYPELIST_H
#define __ITENSOR_TYPELIST_H

//
// Ideas for improvement:
// o Build some tools which do static_asserts
//   on TypeLists e.g. assert that a type is
//   or is not contained
// o Use std::false_type (is there a std::true_type?)
//   instead of NoneType?
// o Write general "find_if" like class
//

namespace itensor {

template<typename... Ts>
struct TypeList { };

template<typename T, typename... Ts>
struct TypeList<T,Ts...> : TypeList<Ts...>
    {
    using Type = T;
    using Next = TypeList<Ts...>;
    };

template<typename TList>
using frontType = typename TList::Type;

template<typename TList>
using popFront = typename TList::Next;

struct NoneType { };

template<typename T, typename ElseType>
using ifTypeElse = std::conditional_t<not std::is_same<T,NoneType>::value,
                                      T,
                                      ElseType>;



template<typename T, typename TList>
struct CheckContainsType : CheckContainsType<T,popFront<TList>>
    {
    using Test = std::conditional_t<std::is_same<frontType<TList>,T>::value,T,NoneType>;
    using ParentResult = typename CheckContainsType<T,popFront<TList>>::Result;
    using Result = ifTypeElse<Test,ParentResult>;
    };
template<typename T>
struct CheckContainsType<T,TypeList<>>
    {
    using Result = NoneType;
    };

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

} //namespace itensor

#endif
