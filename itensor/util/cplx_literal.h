//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_CPLX_LITERAL_H
#define __ITENSOR_CPLX_LITERAL_H

#include "itensor/types.h"

namespace itensor {

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
Cplx inline
operator*(T i, const Cplx& z) { return Real(i)*z; }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
Cplx
operator*(const Cplx& z, T i) { return z*Real(i); }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
Cplx inline
operator/(T i, const Cplx& z) { return Real(i)/z; }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
Cplx
operator/(const Cplx& z, T i) { return z/Real(i); }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
Cplx 
operator+(T i, const Cplx& z) { return Real(i)+z; }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
Cplx 
operator+(const Cplx& z, T i) { return z+Real(i); }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
Cplx 
operator-(T i, const Cplx& z) { return Real(i)-z; }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
Cplx 
operator-(const Cplx& z, int i) { return z-Real(i); }

Cplx constexpr inline
operator"" _i(unsigned long long int i)
    {
    return Cplx(0.,i);
    }

Cplx constexpr inline
operator"" _i(long double x)
    {
    return Cplx(0.,x);
    }

} //namespace itensor


#endif
