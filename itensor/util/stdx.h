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

//
//Useful for disabling candidate template functions
//where if Expression fails to compile (fails template
//substitution), then a lower-precedence overload
//will be selected.
//
template<typename Expression, typename ReturnValue>
using if_compiles_return = ReturnValue;

//
//Dummy argument types to simplify
//template overload precedence.
//Type "select_overload" auto converts to
//choice<1>, then choice<2>, etc.
//in decreasing order of precendence.
//(credit to R. Martinho Fernandes)
//
//Usage:
// funcImpl(..., choice<2>) { }
// funcImpl(..., choice<1>) { }
// func(...) { funcImpl(...,select_overload{}); }
//
template<unsigned I>
struct choice : choice<I+1> { constexpr choice(){} };

template<>
struct choice<10> { constexpr choice(){} };

struct select_overload : choice<1> { constexpr select_overload(){} };


template<typename... VArgs>
auto
make_array(VArgs&&... vargs)
    -> std::array<typename std::common_type<VArgs...>::type,sizeof...(VArgs)>
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
void
fill(Container&& C,
     T&& val)
    {
    std::fill(C.begin(),C.end(),std::forward<T>(val));
    }

template<typename Container,
         typename T>
auto
find(Container&& C,
     T&& val) -> decltype(C.begin())
    {
    return std::find(C.begin(),C.end(),std::forward<T>(val));
    }

template<typename Container,
         class UnaryCmpFunc>
auto
find_if(Container&& C,
        UnaryCmpFunc&& f) -> decltype(C.begin())
    {
    return std::find_if(C.begin(),C.end(),std::forward<UnaryCmpFunc>(f));
    }

} //namespace stdx

#endif
