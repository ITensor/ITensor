//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_COMPLEX_H
#define __ITENSOR_COMPLEX_H

#include "types.h"
#include <complex>

namespace itensor {

using Complex = std::complex<Real>;

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
Complex inline
operator*(T i, const Complex& z) { return Real(i)*z; }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
Complex
operator*(const Complex& z, T i) { return z*Real(i); }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
Complex inline
operator/(T i, const Complex& z) { return Real(i)/z; }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
Complex
operator/(const Complex& z, T i) { return z/Real(i); }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
Complex 
operator+(T i, const Complex& z) { return Real(i)+z; }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
Complex 
operator+(const Complex& z, T i) { return z+Real(i); }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
Complex 
operator-(T i, const Complex& z) { return Real(i)-z; }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
Complex 
operator-(const Complex& z, int i) { return z-Real(i); }

Complex constexpr inline
operator"" _i(unsigned long long int i)
    {
    return Complex(0.,i);
    }

Complex constexpr inline
operator"" _i(long double x)
    {
    return Complex(0.,x);
    }

}; //namespace itensor


#endif
