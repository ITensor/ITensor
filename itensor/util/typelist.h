//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TYPELIST_H
#define __ITENSOR_TYPELIST_H

namespace itensor {

template<typename... Ts>
struct TypeList { };

template<typename T, typename... Ts>
struct TypeList<T,Ts...> : TypeList<Ts...>
    {
    using Type = T;
    using Next = TypeList<Ts...>;
    };


} //namespace itensor

#endif
