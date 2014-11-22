//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_FUNCTIONS_H
#define __ITENSOR_FUNCTIONS_H

namespace itensor {
namespace detail {

template <typename Set1,
          typename Set2,
          typename Perm>
void
calc_permutation(const Set1& s1,
                 const Set2& s2,
                 Perm& p)
    {
    for(size_t i2 = 0; i2 < s2.size(); ++i2)
        {
        const auto& v2 = s2[i2];
        bool found = false;
        for(size_t i1 = 0; i1 < s1.size(); ++i1)
            {
            if(v2 == s1[i1])
                {
                p[i1] = i2;
                found = true;
                break;
                }
            }

        if(!found)
            {
            throw ITError("sets are not permutations of each other");
            }
        }
    }

template <typename Set1,
          typename Set2,
          typename RType,
          typename Map>
void
permute_map(const Set1& s1,
            const Set2& s2,
            RType& r,
            Map&& m)
    {
    for(size_t i2 = 0; i2 < s2.size(); ++i2)
        {
        const auto& v2 = s2[i2];
        bool found = false;
        for(size_t i1 = 0; i1 < s1.size(); ++i1)
            {
            if(v2 == s1[i1])
                {
                r[i1] = m(s2[i2]);
                found = true;
                break;
                }
            }

        if(!found)
            {
            throw ITError("sets are not permutations of each other");
            }
        }
    }

template <typename Container, typename Item>
bool
contains(const Container& C,
         const Item& I)
    {
    for(const auto& c : C) 
        {
        if(I == c) return true;
        }
    return false;
    }


template <typename Ret, class T, typename V>
auto 
call_impl(T&& obj, V&& v, int) -> decltype(obj(v))
    {
    return obj(v);
    }
template <typename Ret, class T, typename V>
Ret
call_impl(T&& obj, V&& v, long) 
    {
    throw std::runtime_error("Object does not support operator() for specified type.");
    return Ret();
    }
//
// The call(obj,v) function uses substitution-failure-is-not-an-error (sfinae)
// to either "plug" v into obj's operator() method if it has one defined
// or else to throw a runtime error if not.
// Use call(obj,v) to convert the absence of a specific operator() method
// to be a run-time error instead of a compile-time error.
//
template <typename Ret, class T, typename V>
Ret
call(T&& obj, V&& v)
    {
    return call_impl<Ret,T,V>(std::forward<decltype(obj)>(obj),v,0);
    }

/////////////////////////

template <class T, typename V>
void 
call_impl(T&& obj, V&& v, int)
    {
    obj(v);
    }
template <class T, typename V>
void
call_impl(T&& obj, V&& v, long) 
    {
    throw std::runtime_error("Object does not support operator() for specified type.");
    }
//
// The call(obj,v) function uses substitution-failure-is-not-an-error (sfinae)
// to either "plug" v into obj's operator() method if it has one defined
// or else to throw a runtime error if not.
// Use call(obj,v) to convert the absence of a specific operator() method
// to be a run-time error instead of a compile-time error.
//
template <class T, typename V>
void
call(T&& obj, V&& v)
    {
    call_impl<T,V>(std::forward<decltype(obj)>(obj),v,0);
    }

/////////////////////////

template <typename Ret, class T, typename V1, typename V2>
auto 
call_impl(T&& obj, V1&& v1, V2&& v2, int) -> decltype(obj(v1,v2))
    {
    return obj(v1,v2);
    }
template <typename Ret, class T, typename V1, typename V2>
Ret
call_impl(T&& obj, V1&& v1, V2&& v2, long) 
    {
    throw std::runtime_error("Object does not support operator() for specified type.");
    return Ret();
    }
//
// The call(obj,v1,v2) function uses substitution-failure-is-not-an-error (sfinae)
// to either "plug" v1,v2 into obj's operator() method if it has one defined
// or else to throw a runtime error if not.
// Use call(obj,v1,v2) to convert the absence of a specific operator() method
// to be a run-time error instead of a compile-time error.
//
template <typename Ret, class T, typename V1, typename V2>
Ret
call(T&& obj, V1&& v1, V2&& v2)
    {
    return call_impl<Ret,T,V1,V2>(std::forward<T>(obj),
                                  std::forward<V1>(v1),
                                  std::forward<V2>(v2),0);
    }

/////////////////////////

template <class T, typename V1, typename V2>
void 
call_impl(T&& obj, V1&& v1, V2&& v2, int)
    {
    obj(v1,v2);
    }
template <class T, typename V1, typename V2>
void
call_impl(T&& obj, V1&& v1, V2&& v2, long) 
    {
    throw std::runtime_error("Object does not support operator() for specified type.");
    }
//
// The call(obj,v1,v2) function uses substitution-failure-is-not-an-error (sfinae)
// to either "plug" v1,v2 into obj's operator() method if it has one defined
// or else to throw a runtime error if not.
// Use call(obj,v1,v2) to convert the absence of a specific operator() method
// to be a run-time error instead of a compile-time error.
//
template <class T, typename V1, typename V2>
void
call(T&& obj, V1&& v1, V2&& v2)
    {
    call_impl<T,V1,V2>(std::forward<T>(obj),
                       std::forward<V1>(v1),
                       std::forward<V2>(v2),0);
    }


}; //namespace detail
}; //namespace itensor

#endif

