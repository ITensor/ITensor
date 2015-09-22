//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_CPLX_LITERAL_H
#define __ITENSOR_CPLX_LITERAL_H

#include <complex>

namespace itensor {


template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
std::complex<double> inline
operator*(T i, std::complex<double> const& z) { return double(i)*z; }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
std::complex<double>
operator*(std::complex<double> const& z, T i) { return z*double(i); }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
std::complex<double> inline
operator/(T i, std::complex<double> const& z) { return double(i)/z; }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
std::complex<double>
operator/(std::complex<double> const& z, T i) { return z/double(i); }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
std::complex<double> 
operator+(T i, std::complex<double> const& z) { return double(i)+z; }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
std::complex<double> 
operator+(std::complex<double> const& z, T i) { return z+double(i); }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
std::complex<double> 
operator-(T i, std::complex<double> const& z) { return double(i)-z; }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
std::complex<double> 
operator-(std::complex<double> const& z, T i) { return z-double(i); }

std::complex<double> constexpr inline
operator"" _i(unsigned long long int i)
    {
    return std::complex<double>(0.,i);
    }

std::complex<double> constexpr inline
operator"" _i(long double x)
    {
    return std::complex<double>(0.,x);
    }

} //namespace itensor


#endif
