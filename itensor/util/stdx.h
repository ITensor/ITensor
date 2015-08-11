//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_STDX
#define __ITENSOR_STDX

#include <array>
#include <vector>
#include <algorithm>

namespace stdx {

template<typename Expression, typename ReturnValue>
using if_compiles_return = ReturnValue;

template<typename... VArgs>
auto
make_array(VArgs&&... vargs)
    -> std::array<std::common_type_t<VArgs...>,sizeof...(VArgs)>
    {
    return {{ std::forward<VArgs>(vargs)... }};
    }

template<typename T>
std::vector<T>
reserve_vector(typename std::vector<T>::size_type size)
    {
    std::vector<T> v;
    v.reserve(size);
    return v;
    }

template<typename Container,
         typename T>
auto
find(Container&& C,
     T&& val)
    {
    return std::find(C.begin(),C.end(),std::forward<T>(val));
    }

template<typename Container,
         class UnaryCmpFunc>
auto
find_if(Container&& C,
        UnaryCmpFunc&& f)
    {
    return std::find_if(C.begin(),C.end(),std::forward<UnaryCmpFunc>(f));
    }

} //namespace stdx

#endif
