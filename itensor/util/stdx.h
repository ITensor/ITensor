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

template<typename T>
using remove_const_t = typename std::remove_const<T>::type;

template<typename T>
using remove_pointer_t = typename std::remove_pointer<T>::type;

template<typename T>
using remove_reference_t = typename std::remove_reference<T>::type;

template<typename T>
using add_pointer_t = typename std::add_pointer<T>::type;

template<typename T>
using decay_t = typename std::decay<T>::type;

template<bool B, class T = void>
using enable_if_t = typename std::enable_if<B,T>::type;

template<bool B, typename T1, typename T2>
using conditional_t = typename std::conditional<B,T1,T2>::type;

template<bool B, typename T1, typename T2>
using conditional_t = typename std::conditional<B,T1,T2>::type;

template<class T>
using result_of_t = typename std::result_of<T>::type;

//
//Useful for disabling candidate template functions
//where if Expressions.. fail to compile (fail template
//substitution), then a lower-precedence overload
//will be selected.
//
template<typename ReturnValue, typename... Expressions>
using if_compiles_return = ReturnValue;

//Helper type for making static_assert always fail,
//but only if a given template is instantiated
//(unlike std::false_type, this type depends on
//the typename T so will not be evaluated until
//the template is instantiated)
template<typename T, typename... Rest>
struct false_regardless_of : public std::false_type
    {
    using ignored_type = T;
    };

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

template<typename Container,
         class UnaryFunc>
void
for_each(Container&& C,
         UnaryFunc&& f)
    {
    std::for_each(std::begin(C),std::end(C),std::forward<UnaryFunc>(f));
    }

template<typename Container,
         class CmpFun>
void
sort(Container && C,
     CmpFun && f)
    {
    std::sort(std::begin(C),std::end(C),std::forward<CmpFun>(f));
    }

} //namespace stdx

#endif
